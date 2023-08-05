#!/usr/bin/env python
import os
import json
import hashlib
import requests
import subprocess
import tarfile
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from modules.common import set_threads, set_out_dir, open_gz
from modules.taxonparse import list_subtree


def validate_md5(file, md5):
    """
    Validate downloaded file using md5
    file: path to file -> Path
    md5: path to md5 checksum -> str
    """
    with open(file, "rb") as f:
        m = hashlib.md5()
        for chunk in iter(lambda: f.read(4096), b""):
            m.update(chunk)
        file_md5 = m.hexdigest()
    return file_md5 == md5


def download_file(url, out_dir, rewrite=False):
    """
    Download file from url to out_dir
    url: url to file -> str
    out_dir: path to output directory -> Path
    return: (returncode, out_file)
    """
    if url.startswith("https://"):
        method = "https"
    elif url.startswith("rsync://"):
        method = "rsync"
    else:
        raise Exception(f"Unknown download method: {url}")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / Path(url).name
    if out_file.is_file() and not rewrite:
        print(f"{out_file} already exists")
        returncode = 0
    else:
        try:
            if method == "https":
                with requests.get(url, stream=True) as r, open(out_file, "wb") as f:
                    r.raise_for_status()
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                returncode = 0
            elif method == "rsync":
                returncode = os.system(
                    f"rsync --copy-links --times --quiet {url} {out_dir}"
                )
            else:
                returncode = 1

        except Exception as e:
            print(f"Failed to download {url}")
            print(e)
            returncode = 1

    return (returncode, out_file)


@set_out_dir
def build_acc2taxid(acc2taxid=None, out_dir=None):
    def accession_to_taxid():
        c = 0
        with open_gz(acc2taxid) as gz_f:
            for line in gz_f:
                c += 1
                line = line.strip(" \n").split("\t")
                if c == 1:
                    continue
                yield (line[0], (int(line[2]),))

    with tempfile.TemporaryDirectory(prefix="acc2taxid_", dir=out_dir) as tmp_dir:
        if acc2taxid:
            acc2taxid = Path(acc2taxid)
        else:
            print("Download acc2taxid file")
            returncode, acc2taxid = download_file(
                url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                out_dir=tmp_dir,
            )
            if returncode:
                raise Exception("Failed to download acc2taxid file")
        trie = marisa_trie.RecordTrie("l", accession_to_taxid())
        marisa_file = out_dir / "acc2taxid.marisa"
        trie.save(marisa_file)
    return marisa_file


@set_threads
def convert_blastdb_to_fasta(
    blastdb, basename, taxids=[], min_len=1, compress=False, threads=0
):
    """
    Convert blastdb to fasta file
    blastdb: path to blastdb -> Path
    basename: basename of output fasta file -> Path
    taxids: list of taxids to extract -> list
    """
    taxids = [str(taxid) for taxid in taxids]
    basename = Path(basename)
    out_dir = basename.parent
    fasta_fpath = basename.parent / f"{basename.name}.fa"
    delim = ">>>"
    extract_cmd = [
        "blastdbcmd",
        "-db",
        blastdb,
        "-outfmt",
        f"%a{delim}%T{delim}%s{delim}%l",
    ]
    if taxids:
        taxidlist = out_dir / f"{basename.name}_taxidlist.txt"
        taxidlist.write_text("\n".join(taxids))
        extract_cmd.extend(["-taxidlist", taxidlist, "-target_only"])
    else:
        extract_cmd.extend(["-entry", "all"])

    min_len = 0 if min_len is None else min_len

    with open(fasta_fpath, "w") as f:
        blastdbcmd_proc = subprocess.Popen(extract_cmd, stdout=subprocess.PIPE)

        while True:
            line = blastdbcmd_proc.stdout.readline().decode()
            if not line:
                break

            seqid, taxid, seq, seqlen = line.split(delim)

            if int(seqlen) >= min_len:
                f.write(f">{seqid}|{taxid}\n{seq}\n")

    if compress:
        compress_proc = subprocess.run(["pigz", "-p", threads, fasta_fpath])
        if compress_proc.returncode:
            raise Exception("Failed to compress fasta file")
        fasta_fpath = fasta_fpath.with_suffix(".fa.gz")

    return fasta_fpath


@set_out_dir
def build_blastdb(taxdump_dir, db_type="nt", out_dir=None, resume=False):
    def batch_get_blastdb(url, out_dir):
        tar_fpath = out_dir / Path(url).name
        md5_fpath = tar_fpath.with_suffix(".gz.md5")

        rewrite = False
        if tar_fpath.is_file() and md5_fpath.is_file():
            rewrite = not validate_md5(tar_fpath, md5_fpath.read_text().split(" ")[0])
            print(f"{tar_fpath} exists, and it is valid")
        else:
            print("Downloding", url)
            returncode, tar_fpath = download_file(url, out_dir, rewrite=rewrite)
            if returncode:
                return (returncode, url)
            md5_url = url + ".md5"
            returncode, md5_fpath = download_file(md5_url, out_dir)
            if returncode:
                return (returncode, url)
            md5 = md5_fpath.read_text().split(" ")[0]
            print("Validating", tar_fpath)
            if not validate_md5(tar_fpath, md5):
                print("MD5 validation failed")
                return (1, url)
            print("Extracting", tar_fpath)
            with tarfile.open(tar_fpath, "r:gz") as tar:
                tar.extractall(out_dir)
        print("Done", tar_fpath)
        return (0, url)

    human_fa = out_dir / f"human_{db_type}.fna"
    non_human_fa = out_dir / f"non_human_{db_type}.fna"
    if human_fa.is_file() and non_human_fa.is_file():
        print("Human and non-human fasta files already exist")
        return {
            "fasta": {"human": human_fa, "non_human": non_human_fa},
        }

    seq_type = "nucl" if db_type == "nt" else "prot"
    metadata_json_url = (
        f"https://ftp.ncbi.nlm.nih.gov/blast/db/{db_type}-{seq_type}-metadata.json"
    )
    returncode, metadata_json = download_file(url=metadata_json_url, out_dir=out_dir)
    assert returncode == 0, "Failed to download nt metadata"

    metadata = json.loads(metadata_json.read_text())
    urls = [url.replace("ftp://", "https://") for url in metadata["files"]]
    blastdb_dir = out_dir / "blastdb"
    args = [(url, blastdb_dir) for url in urls]
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(batch_get_blastdb, *arg) for arg in args]
        results = [future.result() for future in as_completed(futures)]
    if any([result[0] for result in results]):
        raise Exception("Failed to download blastdb")
    db_basename = blastdb_dir / db_type

    all_taxids = list_subtree(taxids=[131567], taxdump_dir=taxdump_dir)
    excluded_taxids = list_subtree(taxids=[33090, 33208], taxdump_dir=taxdump_dir)
    human_taxids = list_subtree(taxids=[9606], taxdump_dir=taxdump_dir)
    non_human_taxids = all_taxids - human_taxids - excluded_taxids

    db_args = [
        (db_basename, non_human_fa.with_suffix(), non_human_taxids, 50, False),
        (db_basename, human_fa.with_suffix(), human_taxids, 50, False),
    ]
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(convert_blastdb_to_fasta, *arg) for arg in db_args]
        results = [future.result() for future in as_completed(futures)]

    if any([result[0] for result in results]):
        raise Exception("Failed to convert blastdb to fasta")

    return {
        "blastdb": db_basename,
        "fasta": {"human": human_fa, "non_human": non_human_fa},
    }


@set_threads
@set_out_dir
def build_human_db(human_fa, out_dir, threads=0):
    """
    Build human db for bowtie2 and HISAT2
    human_nt: path to human nt fasta -> Path
    blastdb_nt: path to nt fasta -> Path
    acc2taxid: path to acc2taxid marisa trie -> Path

    """
    bowtie2_dir = out_dir / "bowtie"
    bowtie_base_index = bowtie2_dir / "nt_human"
    hisat2_url = "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz"
    hisat2_dir = out_dir / Path(hisat2_url).name.replace(".tar.gz", "_hisat")
    hisat2_dir = out_dir / Path(hisat2_url).name.replace(".tar.gz", "_hisat")

    bowtie2_dir.mkdir(parents=True, exist_ok=True)
    if (
        len([f for f in bowtie2_dir.glob("nt_human*")]) > 0
        and len([f for f in hisat2_dir.glob("*.ht2")]) > 0
    ):
        print("Human bowtie2 and hisat2 indexes already exist")
        return {
            "bowtie2": bowtie_base_index,
            "hisat2": hisat2_dir / "genome_tran",
        }

    bowtie_index_p = subprocess.Popen(
        ["bowtie2-build", "--threads", threads, human_fa, bowtie_base_index]
    )

    returncode, hisat2_tar = download_file(url=hisat2_url, out_dir=out_dir)
    if returncode:
        raise Exception("Failed to download hisat2 index")

    tar_dirname = None
    with tarfile.open(hisat2_tar, "r:gz") as tar:
        for tarinfo in tar:
            if tarinfo.isfile() and tarinfo.name.endswith(".ht2"):
                tar.extract(tarinfo, out_dir)
            elif tarinfo.isdir():
                tar_dirname = tarinfo.name
    hisat2_dir = out_dir / tar_dirname
    hisat2_dir = hisat2_dir.rename(out_dir / f"{tar_dirname}_hisat")
    hisat2_tar.unlink()
    bowtie_index_p.wait()

    return {
        "bowtie2": bowtie_base_index,
        "hisat2": hisat2_dir / "genome_tran",
    }


@set_threads
def build_human_hisat2_index(out_dir, threads=0):
    human_fna_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
    human_gtf_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"

    for url in (human_fna_url, human_gtf_url):
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out_dir / Path(url).name, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
    for f in (out_dir / Path(url).name for url in (human_fna_url, human_gtf_url)):
        subprocess.run(["pigz", "-p", threads, f])
    human_fna = (out_dir / Path(human_fna_url).name).with_suffix("")
    human_gtf = (out_dir / Path(human_gtf_url).name).with_suffix("")
    exon_f = out_dir / "exon.ss"
    splice_f = out_dir / "splice.ss"
    hisat_base_index = out_dir / "chm13v2.0_plusY"
    build_proc = subprocess.run(
        [
            "hisat2-build",
            "-p",
            threads,
            "--exon",
            exon_f,
            "--ss",
            splice_f,
            human_fna,
            hisat_base_index,
        ]
    )


@set_threads
@set_out_dir
def build_db(out_dir=None, threads=0):
    taxdump_dir = out_dir / "taxdump"
    taxdump_dir.mkdir(parents=True, exist_ok=True)
    if len([f.is_file() for f in taxdump_dir.glob("*dmp")]) > 0:
        print("Taxdump already exists")
    else:
        returncdoe, taxdump_dir = download_file(
            url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            out_dir=out_dir,
        )
        if returncdoe:
            raise Exception("Failed to download taxdump")
        with tarfile.open(taxdump_dir, "r:gz") as tar:
            tar.extractall(taxdump_dir)
        taxdump_dir.unlink()

    blastdbs = build_blastdb(taxdump_dir=taxdump_dir, out_dir=out_dir, db_type="nt")

    human_dbs = build_human_db(
        human_fa=blastdbs["fasta"]["human"], out_dir=out_dir / "human", threads=0
    )

    return {"blastdb": blastdbs, "human_idx": human_dbs, "taxdump": taxdump_dir}

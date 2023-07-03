#!/usr/bin/env python
import os
import hashlib
import requests
import subprocess
import marisa_trie
import tempfile
import cython
from pathlib import Path

from backend.common import set_threads, set_out_dir, open_gz
from backend.taxonparse import get_lineage


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
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / Path(url).name
    if out_file.is_file() and not rewrite:
        print(f"{out_file} already exists")
        retruncode = 0
    else:
        try:
            with requests.get(url, stream=True) as r, open(out_file, "wb") as f:
                r.raise_for_status()
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            retruncode = 0
        except Exception as e:
            print(f"Failed to download {url}")
            print(e)
            retruncode = 1

    return (retruncode, out_file)


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


@set_out_dir
def build_nt_db(out_dir=None, rewrite=False):
    fasta_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
    print("Start downloading nt fasta file")
    returncode, fasta = download_file(url=fasta_url, out_dir=out_dir, rewrite=rewrite)
    if returncode:
        raise Exception("Failed to download nt fasta file")

    md5_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz.md5"
    returncode, md5_file = download_file(url=md5_url, out_dir=out_dir)
    if returncode:
        raise Exception("Failed to download nt md5 file")
    print("Start validating nt fasta file")
    if validate_md5(file=fasta, md5=md5_file.read_text().split(" ")[0]):
        print("nt fasta file is valid")
        return fasta
    else:
        raise Exception("nt fasta file is invalid")


@set_threads
def convert_blastdb_to_fasta(
    blastdb, basename, taxids=[], min_len=1, compress=True, threads=0
):
    """
    Convert blastdb to fasta file
    blastdb: path to blastdb -> Path
    basename: basename of output fasta file -> Path
    taxids: list of taxids to extract -> list
    """
    taxids = [str(taxid) for taxid in taxids]
    basename = Path(basename)
    fasta_file = basename.parent / f"{basename.name}.fna"
    split_str = "----"
    extract_cmd = [
        "blastdbcmd",
        "-db",
        blastdb,
        "-outfmt",
        f"%a{split_str}%T{split_str}%s",
    ]
    taxid_check = dict([(taxid, True) for taxid in taxids])
    if taxids:
        extract_cmd.extend(["-taxids", ",".join(taxids)])
    else:
        extract_cmd.extend(["-entry", "all"])

    min_len = 0 if min_len is None else min_len

    with open(fasta_file, "w") as f:
        blastdbcmd_proc = subprocess.Popen(extract_cmd, stdout=subprocess.PIPE)

        while True:
            line = blastdbcmd_proc.stdout.readline().decode()
            if not line:
                break

            seqid, taxid, seq = line.strip(" \n").split(split_str)
            if not taxid_check.get(taxid, False):
                continue
            if len(seq) <= min_len:
                continue

            f.write(f">{seqid}|{taxid}\n{seq}\n")

    if compress:
        compress_proc = subprocess.run(["pigz", "-p", threads, fasta_file])
        if compress_proc.returncode:
            raise Exception("Failed to compress fasta file")
        fasta_file = fasta_file.with_suffix(".fna.gz")

    return fasta_file


@set_out_dir
def split_nt(fasta, acc2taxid_marisa, taxdump_dir, out_dir=None, min_len=50):
    """
    Split nt fasta to human and non-human
    fasta: path to nt fasta -> Path
    out_dir: path to output directory -> Path
    """
    stat = {
        "total_n": 0,
        "human_n": 0,
        "non_human_n": 0,
        "removed": {
            "no_taxid": {"n": 0, "seqids": []},
            "less_than": {"cutoff": min_len, "n": 0, "seqids": []},
            "black_taxa": {"n": 0, "seqids": []},
        },
    }

    acc2taxid_trie = Acc2Taxid(acc2taxid_marisa)

    taxids = set()
    c = 0
    removed = []
    positions = {}
    last_pos = 0
    with open_gz(fasta) as f:
        while True:
            line = f.readline()
            if not line:
                break
            offset = f.tell()
            line_len = offset - last_pos
            if line[0] == ">":
                seqid = line.strip(" >\n").split(".")[0]
                taxid = acc2taxid_trie.search_taxid(seqid)
                taxid = str(taxid) if taxid is not None else taxid
                if taxid:
                    taxids.add(taxid)
                else:
                    removed.append(seqid)
                    seqid = None
                positions[seqid] = {
                    "header": {"fr": last_pos, "len": line_len},
                    "seq": {"fr": last_pos + line_len + 1, "len": 0, "seqlen": 0},
                }
                c += 1
                if c % 1000000 == 0:
                    print(f"Processed {c} sequences")
            else:
                if seqid is None:
                    continue
                seqid_end = (
                    positions[seqid]["header"]["fr"] + positions[seqid]["header"]["len"]
                )
                positions[seqid]["seq"]["len"] = offset - seqid_end
                positions[seqid]["seq"]["seqlen"] += line_len - 1

            last_pos = offset

    stat["total_n"] = c
    print(f"Total number of sequences: {c}")

    taxid2lineage = get_lineage(taxids, taxdump_dir)

    conditions = {}
    for taxid, lineage in taxid2lineage.items():
        if lineage.get("species", {}).get("taxid", "") == "9606":
            conditions[taxid] = 0
        elif lineage.get("kingdom", {}).get("taxid", "") in ("33090", "33208"):
            conditions[taxid] = 1
        else:
            conditions[taxid] = 2

    human_fasta = out_dir / "human.fa"
    nhuman_fasta = out_dir / "non_human.fa"
    c = 0
    con = 1
    with open_gz(fasta, "rt") as f, open(human_fasta, "w") as human_f, open(
        nhuman_fasta, "w"
    ) as nhuman_f:
        for seqid, offset_rec in positions.items():
            c += 1
            taxid = offset_rec["taxid"]
            if taxid is None:
                stat["removed"]["no_taxid"]["n"] += 1
                stat["removed"]["no_taxid"]["seqids"].append(seqid)
                continue
            elif offset_rec["seq"]["seqlen"] < min_len:
                stat["removed"]["less_than"]["n"] += 1
                stat["removed"]["less_than"]["seqids"].append(seqid)
                continue

            con = conditions.get(taxid, 2)
            f.seek(offset_rec["header"]["fr"])
            rec_len = offset_rec["header"]["len"] + offset_rec["seq"]["len"]
            if con == 0:
                human_f.write(f.read(rec_len))
                stat["human_n"] += 1
            elif con == 1:
                stat["removed"]["black_taxa"]["seqids"].append(seqid)
            elif con == 2:
                nhuman_f.write(f.read(rec_len))
                stat["non_human_n"] += 1

            if c % 1000000 == 0:
                break

    return stat


@set_threads
def build_human_db(human_nt, out_dir, threads=0):
    """
    Build human db for bowtie2 and HISAT2
    human_nt: path to human nt fasta -> Path
    blastdb_nt: path to nt fasta -> Path
    acc2taxid: path to acc2taxid marisa trie -> Path

    """
    bowtie2_dir = out_dir / "bowtie2"
    hisat2_dir = out_dir / "hisat2"
    bowtie2_dir.mkdir(parents=True, exist_ok=True)
    hisat2_dir.mkdir(parents=True, exist_ok=True)

    # nt_human_fasta = convert_blastdb_to_fasta(
    #     blastdb=blastdb_nt,
    #     basename=bowtie2_dir / "nt_human",
    #     threads=threads,
    #     compress=False,
    #     taxids=["9606"],
    # )
    bowtie_base_index = bowtie2_dir / "human_nt"
    bowtie_index_p = subprocess.Popen(
        ["bowtie2-build", "--threads", threads, human_nt, bowtie_base_index]
    )

    hisat2_url = "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz"
    returncode, hisat2_fname = download_file(url=hisat2_url, out_dir=out_dir)
    if returncode:
        raise Exception("Failed to download hisat2 index")
    tar_proc = subprocess.run(["tar", "-xf", hisat2_fname, "-C", hisat2_dir])
    bowtie_index_p.wait()
    if tar_proc.returncode:
        raise Exception("Failed to extract hisat2 index")

    return {
        "bowtie2": bowtie_base_index,
        "hisat2": hisat2_fname.name.replace(".tar.gz", "") + "genome_tran",
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


class Acc2Taxid:
    def __init__(self, acc2taxid):
        """
        acc2taxid: path to acc2taxid marisa trie -> Path
        """
        self.trie = marisa_trie.RecordTrie("L").mmap(str(acc2taxid))

    def search_taxid(self, acc):
        match = self.trie.get(acc)
        if match:
            return match[0][0]

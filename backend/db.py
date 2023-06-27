#!/usr/bin/env python
import os
import requests
import subprocess
import marisa_trie
import tempfile
from pathlib import Path

from common import set_threads, set_out_dir, open_gz


def download_file(url, out_dir):
    """
    Download file from url to out_dir
    url: url to file -> str
    out_dir: path to output directory -> Path
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / Path(url).name
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(out_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return out_file


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


@set_threads
def build_human_db(blastdb_nt, out_dir, threads=0):
    """
    Build human db for bowtie2 and HISAT2
    """
    bowtie2_dir = out_dir / "bowtie2"
    hisat2_dir = out_dir / "hisat2"
    bowtie2_dir.mkdir(parents=True, exist_ok=True)
    hisat2_dir.mkdir(parents=True, exist_ok=True)

    nt_human_fasta = convert_blastdb_to_fasta(
        blastdb=blastdb_nt,
        basename=bowtie2_dir / "nt_human",
        threads=threads,
        compress=False,
        taxids=["9606"],
    )
    bowtie_base_index = bowtie2_dir / "nt_human"
    bowtie_index_p = subprocess.Popen(
        ["bowtie2-build", "--threads", threads, nt_human_fasta, bowtie_base_index]
    )

    hisat2_url = "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz"
    with requests.get(hisat2_url, stream=True) as r:
        r.raise_for_status()
        with open(hisat2_dir / Path(hisat2_url).name, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    tar_proc = subprocess.run(
        ["tar", "-xf", hisat2_dir / Path(hisat2_url).name, "-C", hisat2_dir]
    )
    bowtie_index_p.wait()

    return {
        "bowtie2": bowtie_base_index,
        "hisat2": hisat2_dir / Path(hisat2_url).name.replace(".tar.gz", "")
        + "genome_tran",
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
    hisat2_index_p = subprocess.Popen(
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
def index_star(reference, out_dir, threads=0):
    index_cmd = [
        "STAR",
        "--runMode",
        "genomeGenerate",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        out_dir,
        "--genomeFastaFiles",
        reference,
    ]
    star_proc = subprocess.run(index_cmd)
    if star_proc.returncode:
        raise Exception("Failed to index reference genomes")
    return out_dir


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
            acc2taxid = download_file(
                url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                out_dir=tmp_dir,
            )
        trie = marisa_trie.RecordTrie("l", accession_to_taxid())
        marisa_file = out_dir / "acc2taxid.marisa"
        trie.save(marisa_file)
    return marisa_file

#!/usr/bin/env python3
import json
import click
from click_option_group import optgroup
import random
import tempfile
import marisa_trie
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from backend.common import set_threads, CONTEXT_SETTINGS, BaseNameType
from backend.taxonparse import get_lineage

# resolve the task that assign reads to taxids


def parse_paf(paf=None, line=None):
    """
    Parse a PAF file into a dictionary of lists.
    paf: path to PAF file -> str or pathlib.Path
    return: dictionary of lists -> dict
    """
    cols = [
        "qname",
        "qlen",
        "qstart",
        "qend",
        "strand",
        "tname",
        "tlen",
        "tstart",
        "tend",
        "nmatch",
        "alnlen",
        "mapq",
    ]

    def _parse_line(line):
        line = line.rstrip("\n").split("\t")

        line_dct = dict([(col, value) for col, value in zip(cols, line[:12])])
        line_dct["qlen"] = int(line_dct["qlen"])
        line_dct["qstart"] = int(line_dct["qstart"])
        line_dct["qend"] = int(line_dct["qend"])
        line_dct["tlen"] = int(line_dct["tlen"])
        line_dct["tstart"] = int(line_dct["tstart"])
        line_dct["tend"] = int(line_dct["tend"])
        line_dct["nmatch"] = int(line_dct["nmatch"])
        line_dct["alnlen"] = int(line_dct["alnlen"])
        line_dct["pident"] = line_dct["nmatch"] / line_dct["alnlen"]
        line_dct["mapq"] = int(line_dct["mapq"])
        for value in line[12:]:
            tag, tag_type, value = value.split(":")
            line_dct[tag] = int(value) if tag_type == "i" else value
        return line_dct

    if line:
        line = _parse_line(line)
        return line
    elif paf:
        with open(paf, "r") as f:
            for line in f.readlines():
                yield _parse_line(line)
    else:
        raise ValueError("Either paf or line must be provided.")


def call_hits(paf=None, accession2taxid_db=None, taxdump_dir=None):
    """
    Assign reads to taxon
    """
    trie = marisa_trie.RecordTrie("L").mmap(accession2taxid_db)
    taxids = set()
    hit_dct = {}
    for aln_rec in parse_paf(paf=paf):
        if aln_rec["tp"] != "P" or aln_rec["alnlen"] < aln_rec["qlen"] * 0.9:
            continue

        read_id = aln_rec["qname"]
        sample_id = read_id.split(".")[0]
        hit_dct.setdefault(sample_id, {})
        accession = aln_rec["tname"].split(".")[0]
        taxid = trie[accession][0][0]
        taxids.add(taxid)
        hit_dct[sample_id].setdefault(read_id, {})
        hit_dct[sample_id][read_id]["alignment"] = aln_rec
        hit_dct[sample_id][read_id]["taxon"] = {"taxid": taxid}

    lineage_dct = get_lineage(taxids=taxids, taxdump_dir=taxdump_dir)
    for sample_id, reads_dct in hit_dct.items():
        for read_id, read_dct in reads_dct.items():
            taxid = read_dct["taxon"]["taxid"]
            read_dct["taxon"]["lineage"] = lineage_dct[str(taxid)]
    return hit_dct


@set_threads
def aln_with_minimap2(
    queries,
    target,
    k=14,
    w=8,
    preset=None,
    work_dir=None,
    threads=0,
    mm2_args={},
):
    """
    Align query to target with minimap2.
    queries: path to query FASTA files -> list or pathlib.Path
    target: path to target FASTA file -> str or pathlib.Path
    k: k-mer size -> int
    w: minimizer window size -> int
    preset: minimap2 preset -> str
    threads: number of threads -> int
    return: path to PAF file -> str or pathlib.Path
    """
    work_dir = Path().cwd() if work_dir is None else Path(work_dir)
    mm2_cmd = ["minimap2"]
    if mm2_args:
        for arg, value in mm2_args.items():
            arg = arg.replace("_", "-")
            arg = f"-{arg}" if len(arg) == 1 else f"--{arg}"
            if isinstance(value, bool):
                if value:
                    mm2_cmd.append(arg)
            else:
                mm2_cmd.extend([arg, value])
    if "t" not in mm2_args:
        mm2_cmd.extend(["-t", threads])
    if "I" not in mm2_args:
        mm2_cmd.extend(["-I", "8G"])
    if preset:
        mm2_cmd.extend(["-x", preset])
    else:
        mm2_cmd.extend(
            [
                "-k",
                str(k),
                "-w",
                str(w),
            ]
        )

    with tempfile.TemporaryDirectory(prefix="mm2_", dir=work_dir) as tmp_dir:
        mm2_cmd.extend(["--split-prefix", f"{tmp_dir}/mm2"])
        out_paf = f"{work_dir}/mm2.paf"
        mm2_cmd.extend(["-o", out])
        mm2_cmd.append(target)
        if isinstance(queries, list):
            mm2_cmd.extend(queries)
        else:
            mm2_cmd.append(queries)
        print(" ".join([str(i) for i in mm2_cmd]))
        mm2_proc = subprocess.run(mm2_cmd)
        if mm2_proc.returncode != 0:
            raise subprocess.CalledProcessError(mm2_proc.returncode, mm2_proc.args)

    return out_paf


def select_best_alignment(alignments=[], paf=None):
    """
    Select the best alignment from a list of alignments.
    iter_alignments: list of alignments -> list
    return: best alignment -> dict
    """
    best_dct = {}

    for alignment in alignments if alignments else parse_paf(paf=paf):
        if alignment["alnlen"] < alignment["qlen"] * 0.9:
            continue
        qname = alignment["qname"]
        tname = alignment["tname"]
        pident = alignment["pident"]
        if qname not in best_dct or pident > best_dct[qname]["pident"]:
            best_dct[qname] = {"pident": pident, "tnames": [tname], "aln": alignment}
        elif pident == best_dct[qname]["pident"]:
            best_dct[qname]["tnames"].append(tname)

    return best_dct


def reassign_alignments(best_dct):
    """
    Reassign alignments to the best reference.
    best_dct: dictionary of best alignments -> dict
    return: dictionary of re-assigned alignments -> dict
    """
    reassign_dct = {}
    target_dct = {}
    for qname, match_dct in best_dct.items():
        for tname in match_dct["tnames"]:
            target_dct.setdefault(tname, 0)
            target_dct[tname] += 1

    for qname, match_dct in best_dct.items():
        if len(match_dct["tnames"]) == 1:
            reassign_dct[qname] = {
                "best": match_dct["tnames"][0],
                "pident": match_dct["pident"],
                "aln": match_dct["aln"],
            }
        else:
            random.shuffle(match_dct["tnames"])
            reassign_dct[qname] = {
                "best": max(match_dct["tnames"], key=lambda x: target_dct[x]),
                "pident": match_dct["pident"],
                "aln": match_dct["aln"],
            }

    return reassign_dct


def main(queries, reference, out_dir, threads, read_type, paf=None):
    preset = None if read_type == "illumina" else "map-ont"
    out_dir.mkdir(parents=True, exist_ok=True)
    paf = (
        aln_with_minimap2(
            target=reference,
            queries=queries,
            preset=preset,
            work_dir=out_dir,
            threads=threads,
        )
        if paf is None
        else paf
    )
    best_aln_dct = select_best_alignment(paf=paf)
    reassign_dct = reassign_alignments(best_aln_dct)
    aln_json = out_dir / "alignment.json"

    aln_json.write_text(json.dumps(reassign_dct, indent=4))
    return aln_json


@click.command(
    help="Reads mapping to reference genomes", context_settings=CONTEXT_SETTINGS
)
@optgroup.group("Input options")
@optgroup.option(
    "--queries",
    "-q",
    nargs="+",
    help="reads (multiple files allowed)",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
)
@optgroup.option(
    "--reference",
    "-r",
    help="reference genomes",
    type=BaseNameType(),
)
@optgroup.option(
    "--paf",
    "-p",
    help="paf file (mutually exclusive with queries)",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
)
@optgroup.option(
    "--threads",
    "-t",
    help="number of threads",
    type=click.INT,
    default=0,
    show_default=True,
)
@optgroup.option(
    "--read_type",
    "-rt",
    help="read type",
    type=click.Choice(["illumina", "ont"]),
    default="illumina",
    show_default=True,
)
@optgroup.group("Output options")
@optgroup.option(
    "--out_dir",
    "-o",
    help="output directory",
    type=click.Path(exists=False, dir_okay=True, resolve_path=True, path_type=Path),
    required=True,
)
def cli(queries, reference, paf, out_dir, threads, read_type):
    if queries is None and paf is None:
        raise ValueError("Either queries or paf must be provided.")
    if queries is not None and paf is not None:
        raise ValueError("Only one of queries or paf must be provided.")
    if queries is None and reference is None and paf is None:
        raise ValueError("reference and queries must be provided.")

    aln_json = main(
        queries=queries,
        reference=reference,
        paf=paf,
        out_dir=out_dir,
        threads=threads,
        read_type=read_type,
    )
    click.echo(f"Alignment json written to {aln_json}")

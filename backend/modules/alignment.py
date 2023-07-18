#!/usr/bin/env python3
import json
import click
from click_option_group import optgroup, GroupedOption
import random
import tempfile
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from modules.common import set_threads, CONTEXT_SETTINGS, BaseNameType
from modules.taxonparse import get_lineage


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
        mm2_cmd.extend(["-o", out_paf])
        mm2_cmd.append(target)
        if isinstance(queries, (list, tuple)):
            mm2_cmd.extend(queries)
        else:
            mm2_cmd.append(queries)
        print(" ".join([str(i) for i in mm2_cmd]))
        mm2_proc = subprocess.run(mm2_cmd)
        if mm2_proc.returncode != 0:
            raise subprocess.CalledProcessError(mm2_proc.returncode, mm2_proc.args)

    return out_paf


def call_hits(paf=None, taxdump_dir=None):
    """
    Assign reads to taxon
    """
    # taxids = set()
    sample2hits = {}
    for aln_rec in parse_paf(paf=paf):
        if aln_rec["tp"] != "P" or aln_rec["alnlen"] < aln_rec["qlen"] * 0.9:
            continue

        read_id = aln_rec["qname"]
        sample_id = read_id.split(".")[0]
        sample2hits.setdefault(sample_id, [])
        sample2hits[sample_id].append(aln_rec)
    return sample2hits


def select_best_alignment(alignments=[], paf=None):
    """(Deprecated)
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
    """(Deprecated)
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
    sample2hits = call_hits(paf=paf)
    hit_jsons = []
    for sample_id, hits in sample2hits.items():
        hit_json = out_dir / f"{sample_id}.hit.json"
        hit_json.write_text(json.dumps(hits, indent=4))
    return hit_jsons


class OptionEatAll(GroupedOption):
    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop("save_other_options", True)
        nargs = kwargs.pop("nargs", -1)
        assert nargs == -1, "nargs, if set, must be -1 not {}".format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


@click.command(
    help="Reads mapping to reference genomes", context_settings=CONTEXT_SETTINGS
)
@optgroup.group("Input options")
@optgroup.option(
    "--queries",
    "-q",
    nargs=-1,
    cls=OptionEatAll,
    help="reads (multiple files allowed)",
    type=tuple,
)
@optgroup.option(
    "--reference",
    "-r",
    help="reference genomes",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
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

    hit_jsons = main(
        queries=queries,
        reference=reference,
        paf=paf,
        out_dir=out_dir,
        threads=threads,
        read_type=read_type,
    )
    for hit_json in hit_jsons:
        click.echo(f"hit json written to {hit_json}")

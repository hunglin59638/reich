#!/usr/bin/env python
import click
import tempfile
import subprocess
from pathlib import Path

from common import (
    set_threads,
    set_out_dir,
    get_basename,
    BaseNameType,
    CONTEXT_SETTINGS,
)


@set_threads
@set_out_dir
def qc_filter(fastq_1, fastq_2=None, threads=0, out_dir=None):
    """
    run fastp to do quality control
    fastq_1: path to fastq file -> Path
    fastq_2: path to fastq file -> Path
    """
    fastq_1_basename = get_basename(fastq_1)
    fastq_2_basename = get_basename(fastq_2) if fastq_2 is not None else None
    qc_cmd = [
        "fastp",
        "-w",
        threads,
        "--dont_eval_duplication",
        "--length_required",
        "35",
        "--qualified_quality_phred",
        "17",
        "--unqualified_percent_limit",
        "15",
        "--n_base_limit",
        "1",
        "--low_complexity_filter",
        "--complexity_threshold",
        "30",
        "--detect_adapter_for_pe",
        "--json",
        out_dir / "fastp.json",
        "--html",
        out_dir / "fastp.html",
        "-i",
        fastq_1,
        "-o",
        out_dir / f"{fastq_1_basename}.qc.fq",
    ]
    if fastq_2 is not None:
        qc_cmd.extend(["-I", fastq_2, "-O", out_dir / f"{fastq_2_basename}.qc.fq"])
    qc_proc = subprocess.run(qc_cmd)
    if qc_proc.returncode:
        raise Exception("Failed to run fastp")
    return {
        "fastq_1": out_dir / f"{fastq_1_basename}.qc.fq",
        "fastq_2": out_dir / f"{fastq_2_basename}qc.fq"
        if fastq_2 is not None
        else None,
        "json": out_dir / "fastp.json",
        "html": out_dir / "fastp.html",
    }


@set_threads
@set_out_dir
def remove_human_reads(
    fastq_1,
    fastq_2=None,
    bowtie2_idx=None,
    hisat2_idx=None,
    threads=0,
    out_dir=None,
):
    """
    Remove human reads from fastq files
    fastq_1: path to fastq file -> Path
    fastq_2: path to fastq file -> Path
    """
    if bowtie2_idx is None and hisat2_idx is None:
        raise Exception("Either bowtie2_idx or hisat2_idx must be provided")

    fastq_1_basename = get_basename(fastq_1)
    fastq_2_basename = get_basename(fastq_2) if fastq_2 is not None else None

    with tempfile.TemporaryDirectory(
        prefix="human_removal", dir=fastq_1.parent
    ) as tmp_dir:
        tmp_dir = Path(tmp_dir)
        for aln_prog, ref_idx in zip(["bowtie2", "hisat2"], [bowtie2_idx, hisat2_idx]):
            if ref_idx is None:
                continue
            aln_cmd = [aln_prog, "-p", threads, "-x", ref_idx]
            if fastq_2 is None:
                aln_cmd.extend(["-U", fastq_1])
            else:
                aln_cmd.extend(["-1", fastq_1, "-2", fastq_2])

            aln_proc = subprocess.Popen(aln_cmd, stdout=subprocess.PIPE)

            if fastq_2 is not None:
                sort_cmd = ["samtools", "sort", "-@", "4", "-n", "-O", "SAM"]
                sort_proc = subprocess.Popen(
                    sort_cmd, stdin=aln_proc.stdout, stdout=subprocess.PIPE
                )

            fastq_cmd = ["samtools", "fastq", "-f", "4", "-@", "4", "-"]
            if fastq_2 is None:
                fastq_1 = tmp_dir / "filter.fq"
                fastq_cmd.extend(["-0", str(fastq_1), "-s", "/dev/null"])
                fastq_proc = subprocess.run(fastq_cmd, stdin=aln_proc.stdout)
                if not fastq_1.is_file():
                    raise Exception("Failed to extract non-human reads")
            else:
                fastq_1 = tmp_dir / "filter_1.fq"
                fastq_2 = tmp_dir / "filter_2.fq"
                fastq_cmd.extend(
                    [
                        "-1",
                        str(fastq_1),
                        "-2",
                        str(fastq_2),
                        "-0",
                        "/dev/null",
                        "-s",
                        "/dev/null",
                        "-n",
                    ]
                )
                fastq_proc = subprocess.run(fastq_cmd, stdin=sort_proc.stdout)

            if fastq_proc.returncode:
                print(" ".join([str(i) for i in fastq_proc.args]))
                raise Exception("Failed to extract non-human reads")
            if fastq_1.is_file():
                fastq_1 = fastq_1.rename(f"{fastq_1_basename}.nonhuman.fq")
            if fastq_2 is not None and fastq_2.is_file():
                fastq_2 = fastq_2.rename(f"{fastq_2_basename}.nonhuman.fq")
    return fastq_1, fastq_2


@set_threads
@set_out_dir
def main(
    fastq_1, fastq_2=None, bowtie2_idx=None, hisat2_idx=None, threads=0, out_dir=None
):
    qc_out = qc_filter(fastq_1, fastq_2, threads=threads, out_dir=out_dir)
    fastq_1 = qc_out["fastq_1"]
    fastq_2 = qc_out["fastq_2"]
    if bowtie2_idx is not None and hisat2_idx is not None:
        nhuman_fastq_1, nhuman_fastq_2 = remove_human_reads(
            fastq_1,
            fastq_2,
            bowtie2_idx=bowtie2_idx,
            hisat2_idx=hisat2_idx,
            threads=threads,
            out_dir=out_dir,
        )
    else:
        raise Exception("bowtie2_idx and hisat2_idx must be provided")
    fastq_1.unlink()
    if fastq_2 is not None:
        fastq_2.unlink()
    return nhuman_fastq_1, nhuman_fastq_2


@click.command(
    help="Remove host reads from fastq files", context_settings=CONTEXT_SETTINGS
)
@click.option(
    "--fastq_1",
    "-1",
    help="read1 fastq file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
)
@click.option(
    "--fastq_2",
    "-2",
    help="read2 fastq file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
)
@click.option(
    "--bowtie2_idx",
    "-b",
    help="bowtie2 index of human",
    type=BaseNameType(),
    required=True,
)
@click.option(
    "--hisat2_idx",
    help="hisat2 index of human",
    type=BaseNameType(),
    required=True,
)
@click.option(
    "--threads",
    "-t",
    help="number of threads, default: 0 (use all available cores)",
    type=int,
    default=0,
    show_default=False,
)
@click.option(
    "--out_dir",
    "-o",
    help="output directory",
    type=click.Path(exists=False, dir_okay=True, resolve_path=True, path_type=Path),
    default=Path().cwd(),
    show_default=True,
)
def cli(fastq_1, fastq_2, bowtie2_idx, hisat2_idx, threads, out_dir):
    main(
        fastq_1=fastq_1,
        fastq_2=fastq_2,
        bowtie2_idx=bowtie2_idx,
        hisat2_idx=hisat2_idx,
        threads=threads,
        out_dir=out_dir,
    )

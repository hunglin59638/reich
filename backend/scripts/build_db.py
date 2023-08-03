#!/usr/bin/env python3
import re
import sys
import click
from pathlib import Path

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
from modules.common import CONTEXT_SETTINGS
from modules.db import build_db


def update_nf_config(db_paths):
    """
    update db path to nextflow.config
    """
    config_fpath = backend_dir / "nextflow_template.config"
    updated_fpath = backend_dir / "nextflow.config"
    bowtie2_idx = db_paths["human_idx"]["bowtie2"]
    hisat2_idx = db_paths["human_idx"]["hisat2"]
    nonhuman_db = db_paths["blastdb"]["fasta"]["non_human"]
    taxdump_dir = db_paths["taxdump"]

    with open(config_fpath, "r") as f_in, open(updated_fpath, "w") as f_out:
        for row in f_in:
            if re.search(r"bowtie2_idx = .+", row):
                row = re.sub(r"bowtie2_idx = .+", f'bowtie2_idx = "{bowtie2_idx}"', row)
            elif re.search(r"hisat2_idx = .+", row):
                row = re.sub(r"hisat2_idx = .+", f'hisat2_idx = "{hisat2_idx}"', row)
            elif re.search(r"nonhuman_db = .+", row):
                row = re.sub(r"nonhuman_db = .+", f'nonhuman_db = "{nonhuman_db}"', row)
            elif re.search(r"taxdump_dir = .+", row):
                row = re.sub(r"taxdump_dir = .+", f'taxdump_dir = "{taxdump_dir}"', row)
            f_out.write(row)
    print(f"updated {updated_fpath}")


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--out_dir",
    "-o",
    help="output directory",
    type=click.Path(exists=False, dir_okay=True, resolve_path=True, path_type=Path),
)
@click.option(
    "--threads",
    "-t",
    help="number of threads, default: 0 (use all available cores)",
    type=int,
    default=0,
)
def main(out_dir, threads):
    db_paths = build_db(out_dir=out_dir, threads=threads)
    update_nf_config(db_paths)


if __name__ == "__main__":
    main()

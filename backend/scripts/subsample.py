#!/usr/bin/env python3
import sys
from pathlib import Path
import subprocess
import click

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modules.common import set_out_dir, CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--reads",
    help="Path to reads file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
)
@click.option(
    "--subsample_n",
    help="Subsample reads to this number",
    type=int,
    required=True,
    default=1000000,
)
@click.option("--out_dir", help="Path to output directory", type=Path, required=True)
@set_out_dir
def main(reads, subsample_n, out_dir):
    i = 0
    line_n = 0
    sample_id = reads.name.split(".")[0]
    print(f"sample id: {sample_id}")
    cmd = ["seqtk", "sample", reads, str(subsample_n)]
    with open(out_dir / f"{sample_id}.subsampled.fq", "w") as f:
        for line in subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout:
            line = line.decode("utf-8")
            line_n += 1
            if line_n % 4 == 1:
                i += 1
                read_id = f"{sample_id}.{i}"
                line = f"@{read_id}\n"

            f.write(line)
    print("output written to", out_dir / f"{sample_id}.subsampled.fq")


if __name__ == "__main__":
    main()

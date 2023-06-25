#!/usr/bin/env python
import sys
import json
import click
from pathlib import Path
from time import sleep
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from backend.common import check_seq_format, set_out_dir, AdvancedJSONEncoder


@click.command(help="Validate the input files")
@click.option(
    "--reads",
    help="input reads",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    required=True,
)
@click.option("--out_dir", "-o", help="output directory")
@set_out_dir
def main(reads, out_dir):
    valid_json = out_dir/"valid.json"
    valid_dct = {"valid_json": valid_json}
    if check_seq_format(reads) != "fastq":
        click.echo(f"{reads} is not a fastq")
        raise click.Abort(f"{reads} is not a fastq")
    valid_dct.update({"valid_fastq": reads})
    valid_json.write_text(json.dumps(valid_dct, indent=4, cls=AdvancedJSONEncoder))
    click.echo("OK")
    sleep(3)


if __name__ == "__main__":
    main()

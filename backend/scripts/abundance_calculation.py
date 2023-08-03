#!/usr/bin/env python3

import sys
import json
import click
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modules.common import set_out_dir, CONTEXT_SETTINGS
from modules.abundance import get_rpm


@click.command(help="Relative abundance calculation", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--taxon_json",
    "-t",
    help="input taxon json",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    required=True,
)
@click.option(
    "--out_dir",
    "-o",
    help="output directory",
    type=click.Path(exists=False, dir_okay=True, resolve_path=True, path_type=Path),
    required=True,
)
@set_out_dir
def main(taxon_json, out_dir):
    read2taxon = json.loads(taxon_json.read_text())

    rpm_dct = get_rpm(read2taxon)
    sample_id = taxon_json.name.split(".")[0]
    out_json = out_dir / f"{sample_id}.rpm.json"
    out_json.write_text(json.dumps(rpm_dct, indent=2))
    click.echo(f"Output: {out_json}")
    return 0


if __name__ == "__main__":
    main()

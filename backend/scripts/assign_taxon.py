#!/usr/bin/env python3
import sys
import json
import click
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modules.common import set_out_dir, CONTEXT_SETTINGS
from modules.taxonparse import assign_taxon


@click.command(help="Assign reads to taxon from paf", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--hit_json",
    help="input {sample_id}.hits.json",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    required=True,
)
@click.option(
    "--taxdump_dir",
    help="taxdump directory",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True, path_type=Path),
    required=True,
)
@click.option("--out_dir", "-o", help="output directory")
@set_out_dir
def main(hit_json, taxdump_dir, out_dir):
    sample_id = hit_json.name.split(".")[0]
    qname2lineage = assign_taxon(hit_json=hit_json, taxdump_dir=taxdump_dir)
    lineage_json = out_dir / f"{sample_id}.taxonomy.json"
    lineage_json.write_text(json.dumps(qname2lineage, indent=4))
    click.echo(f"Output: {lineage_json}")


if __name__ == "__main__":
    main()

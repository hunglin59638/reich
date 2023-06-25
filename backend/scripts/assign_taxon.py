#!/usr/bin/env python3
import sys
import json
import click
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from backend.common import set_out_dir, CONTEXT_SETTINGS
from backend.alignment import call_hits


@click.command(help="Assign reads to taxon from paf", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--paf",
    help="input paf",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    required=True,
)
@click.option(
    "--taxdump_dir",
    help="taxdump directory",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True, path_type=Path),
    required=True,
)
@click.option(
    "--accession2taxid_db",
    "-a",
    help="accession2taxid database of marisa_trie",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True, path_type=Path),
    required=True,
)
@click.option("--out_dir", "-o", help="output directory")
@set_out_dir
def main(paf, taxdump_dir, accession2taxid_db, out_dir):
    hits_dct = call_hits(
        paf=paf, taxdump_dir=taxdump_dir, accession2taxid_db=str(accession2taxid_db)
    )
    for sample_id, reads_dct in hits_dct.items():
        out_json = out_dir / f"{sample_id}.hits.json"
        out_json.write_text(json.dumps(reads_dct, indent=4))
        click.echo(f"Output: {out_json}")


if __name__ == "__main__":
    main()

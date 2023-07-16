#!/usr/bin/env python3
import sys
import json
import click
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from backend.common import set_out_dir, CONTEXT_SETTINGS
from backend.alignment import call_hits
from backend.taxonparse import assign_taxon


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
    # hits_dct = call_hits(
    #     paf=paf, taxdump_dir=taxdump_dir, accession2taxid_db=str(accession2taxid_db)
    # )
    # taxonomy_dct = {}
    # for sample_id, reads_dct in hits_dct.items():
    #     hit_json = out_dir / f"{sample_id}.hits.json"
    #     hit_json.write_text(json.dumps(reads_dct, indent=4))
    #     click.echo(f"Output: {hit_json}")

    #     for read_id, value_dct in reads_dct.items():
    #         taxon_dct = value_dct["taxon"]
    #         aln_dct = value_dct["alignment"]
    #         sp_taxon = taxon_dct["lineage"].get("species")
    #         if not sp_taxon:
    #             continue
    #         sp_taxid = sp_taxon["taxid"]
    #         taxonomy_dct.setdefault(
    #             sp_taxid,
    #             {
    #                 "species_taxon": sp_taxon,
    #                 "genus_taxon": taxon_dct["lineage"].get("genus"),
    #                 "family_taxon": taxon_dct["lineage"].get("family"),
    #                 "alignment_length": [],
    #                 "count": 0,
    #                 "percent_identity": [],
    #             },
    #         )

    #         taxonomy_dct[sp_taxid]["count"] += 1
    #         taxonomy_dct[sp_taxid]["alignment_length"].append(aln_dct["alnlen"])
    #         taxonomy_dct[sp_taxid]["percent_identity"].append(aln_dct["pident"])

    #     for sp_taxid in list(taxonomy_dct.keys()):
    #         taxonomy_dct[sp_taxid]["alignment_length"] = sum(
    #             taxonomy_dct[sp_taxid]["alignment_length"]
    #         ) / len(taxonomy_dct[sp_taxid]["alignment_length"])
    #         taxonomy_dct[sp_taxid]["percent_identity"] = sum(
    #             taxonomy_dct[sp_taxid]["percent_identity"]
    #         ) / len(taxonomy_dct[sp_taxid]["percent_identity"])
    #     taxonomy_json = out_dir / f"{sample_id}.taxonomy.json"
    #     taxonomy_json.write_text(json.dumps(taxonomy_dct, indent=4))
    #     click.echo(f"Output: {taxonomy_json}")


if __name__ == "__main__":
    main()

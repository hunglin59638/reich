#!/usr/bin/env python3
import sys
import json
import click
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modules.common import set_out_dir, CONTEXT_SETTINGS


@click.command(help="Relative abundance calculation", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--hit_json",
    "-i",
    help="input hit json",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    required=True,
)
@click.option(
    "--taxon_json",
    "-t",
    help="input taxonomy json",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    required=True,
)
@click.option(
    "--rpm_json",
    "-r",
    help="input rpm abundance json",
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
def main(hit_json, taxon_json, rpm_json, out_dir):
    cols = ["Taxon", "Score", "Z score", "rPM", "r", "%id", "L"]
    hits = json.loads(hit_json.read_text())
    read2hit = {hit["qname"]: hit for hit in hits}
    read2taxon = json.loads(taxon_json.read_text())
    taxid2sp_taxid = {
        taxon["taxid"]: taxon["lineage"]["species"]["taxid"]
        for read_id, taxon in read2taxon.items()
    }
    taxid2rpms = json.loads(rpm_json.read_text())

    sp_taxid2idents = {}
    sp_taxid2aln_lens = {}
    for read_id, hit in read2hit.items():
        taxid = hit["tname"].split("|")[1]
        sp_taxid = taxid2sp_taxid[taxid]
        ident = hit["pident"]
        aln_len = hit["alnlen"]
        sp_taxid2idents.setdefault(sp_taxid, []).append(ident)
        sp_taxid2aln_lens.setdefault(sp_taxid, []).append(aln_len)

    rows = []

    for sp_taxid, values in taxid2rpms.items():
        rpm = values["rpm"]
        r = values["hit_n"]
        pident = round(
            sum(sp_taxid2idents[sp_taxid]) / len(sp_taxid2idents[sp_taxid]) * 100, 2
        )
        aln_len = round(
            sum(sp_taxid2aln_lens[sp_taxid]) / len(sp_taxid2aln_lens[sp_taxid]), 2
        )
        sp_name = values["taxon"]["lineage"]["species"]["name"]
        rows.append([sp_name, "-", "-", rpm, r, pident, aln_len])

    sample_id = rpm_json.name.split(".")[0]
    report_df = pd.DataFrame(rows, columns=cols)
    report_tsv = out_dir / f"{sample_id}.report.tsv"
    report_df.to_csv(report_tsv, sep="\t", index=False)
    click.echo(f"Output: {report_tsv}")


if __name__ == "__main__":
    main()

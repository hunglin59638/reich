#!/usr/bin/env python3

import os
import subprocess


def get_lineage(taxids, taxdump_dir):
    """
    convert taxid to lineage with taxonkit
    taxids: list of taxids -> list
    taxdump_dir: path to taxdump directory -> Path
    """
    taxids = set([str(taxid) for taxid in taxids])
    lineage_cmd = [
        "taxonkit",
        "lineage",
        "--show-lineage-ranks",
        "--show-lineage-taxids",
        "--data-dir",
        taxdump_dir,
    ]
    lineage_proc = subprocess.run(
        lineage_cmd, input="\n".join(taxids).encode(), stdout=subprocess.PIPE
    )
    if lineage_proc.returncode:
        raise Exception("Failed to convert taxid to lineage")

    lineage_dct = {}
    rank_dct = dict(
        [
            (r, True)
            for r in (
                "superkingdom",
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            )
        ]
    )
    for row in lineage_proc.stdout.decode("utf-8").strip().split("\n"):
        row = row.split("\t")
        if len(row) == 4:
            taxid, names, taxids, ranks = row
            lineage_dct[taxid] = dict(
                [
                    (rank, {"name": name, "taxid": taxid})
                    for name, taxid, rank in zip(
                        names.split(";"), taxids.split(";"), ranks.split(";")
                    )
                    if rank_dct.get(rank)
                ]
            )
        else:
            lineage_dct[row[0]] = []
    return lineage_dct

#!/usr/bin/env python3
import json
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


def list_subtree(taxids, taxdump_dir):
    """
    list all taxids in subtree
    taxids: list of taxids -> list
    taxdump_dir: path to taxdump directory -> Path
    """
    taxids = set([str(taxid) for taxid in taxids])
    list_cmd = [
        "taxonkit",
        "list",
        "--data-dir",
        taxdump_dir,
        "--ids",
        ",".join(taxids),
        "-j",
        "4",
    ]
    list_proc = subprocess.Popen(list_cmd, stdout=subprocess.PIPE)
    sub_taxids = set()
    for row in list_proc.stdout:
        row = row.decode("utf-8").strip()
        if row:
            sub_taxids.add(row)
    return sub_taxids


def assign_taxon(hit_json, taxdump_dir):
    taxids = set()
    hits = json.loads(hit_json.read_text())
    qname2lineage = {}
    for hit in hits:
        qname = hit["qname"]
        acc, taxid = hit["tname"].split("|")
        taxids.add(taxid)
        qname2lineage[qname] = {"taxid": taxid, "lineage": None}

    taxid2lineage = get_lineage(taxids, taxdump_dir)
    for qname in list(qname2lineage.keys()):
        taxid = qname2lineage[qname]["taxid"]
        qname2lineage[qname]["lineage"] = taxid2lineage.get(taxid)
    return qname2lineage

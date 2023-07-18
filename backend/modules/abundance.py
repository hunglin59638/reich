#!/usr/bin/env python


def get_rpm(read2taxon, rank="species"):
    """
    Calculate relative abundance of each taxon

    """
    taxid2rpm = {}
    total_n = 0
    for read_id, taxon in read2taxon.items():
        if taxon["lineage"].get(rank):
            total_n += 1
            taxid = taxon["lineage"][rank]["taxid"]
            taxid2rpm.setdefault(taxid, {"rpm": 0, "hit_n": 0, "taxon": taxon})
            taxid2rpm[taxid]["hit_n"] += 1

    for taxid in list(taxid2rpm.keys()):
        taxid2rpm[taxid]["rpm"] = (taxid2rpm[taxid]["hit_n"] / total_n) * (10**6)
    taxid2rpm = dict(sorted(taxid2rpm.items(), key=lambda x: x[1]["rpm"], reverse=True))
    return taxid2rpm


def get_zscore(taxon, background_model):
    pass

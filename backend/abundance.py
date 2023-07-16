#!/usr/bin/env python


def get_rpm(read2taxon, rank="species"):
    """
    Calculate relative abundance of each taxon

    """
    rpm_dct = {}
    total_n = 0
    for read_id, taxon_dct in read2taxon.items():
        if taxon_dct["lineage"].get(rank):
            total_n += 1
            taxon_taxid = taxon_dct["lineage"][rank]["taxid"]
            rpm_dct.setdefault(taxon_taxid, {"rpm": 0, "hit_n": 0, "taxon": taxon_dct})
            rpm_dct[taxon_taxid]["hit_n"] += 1

    for taxon_taxid in list(rpm_dct.keys()):
        rpm_dct[taxon_taxid]["rpm"] = (rpm_dct[taxon_taxid]["hit_n"] / total_n) * (
            10**6
        )
    rpm_dct = dict(sorted(rpm_dct.items(), key=lambda x: x[1]["rpm"], reverse=True))
    return rpm_dct


def get_zscore(taxon_dct, background_model):
    pass

#!/usr/bin/env python


def get_rpm(read2taxon_dct, rank="species"):
    """
    Calculate relative abundance of each taxon

    """
    rpm_dct = {}
    total_n = 0
    for read_id, taxon_dct in read2taxon_dct.items():
        if taxon_dct["lineage"].get(rank):
            total_n += 1
            taxon_name = taxon_dct["lineage"][rank]["name"]
            rpm_dct.setdefault(taxon_name, {"rpm": 0, "hit_n": 0, "taxon": taxon_dct})
            rpm_dct[taxon_name]["hit_n"] += 1

    for taxon_name in list(rpm_dct.keys()):
        rpm_dct[taxon_name]["rpm"] = (rpm_dct[taxon_name]["hit_n"] / total_n) * (
            10**6
        )
    rpm_dct = dict(sorted(rpm_dct.items(), key=lambda x: x[1]["rpm"], reverse=True))
    print(total_n)
    return rpm_dct


def get_zscore(taxon_dct, background_model):
    pass

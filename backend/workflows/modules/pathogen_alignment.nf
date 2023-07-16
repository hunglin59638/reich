#!/usr/bin/env nextflow 

process pathogen_alignment {
    publishDir "${params.out_dir}/hit", mode: 'copy', pattern: '*.hit.json'

    label "performance"
    input:
        path(reads)
    output:
        path('hit/*.hit.json'), emit: hit_json
    script:
    """
    python $workflow.projectDir/../scripts/pathogen_alignment.py --out_dir hit --queries $reads --reference $params.nonhuman_db --threads $params.threads
    """
}

process assign_taxon {
    publishDir "${params.out_dir}/taxon", mode: 'copy', pattern: '*.taxonomy.json'

    label "normal"
    input:
        path(hit_json)
    output:
        path(hit_json), emit: hit_json
        path("taxon/*.taxonomy.json"), emit: taxon_json
    script:
    """
    python $workflow.projectDir/../scripts/assign_taxon.py --hit_json $hit_json --out_dir taxon --taxdump_dir $params.taxdump_dir
    """
}
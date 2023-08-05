#!/usr/bin/env nextflow 
nextflow.enable.dsl = 2

process validate_reads {
    publishDir "${params.out_dir}" , mode: 'copy', pattern: 'valid_reads/*.json'
    label "normal"
    input:
        path(reads)
    output:
        path('valid_reads/*.json'), emit: valid_json
        path(reads), emit: valid_reads
    script:
    """
    python $workflow.projectDir/scripts/validate_input.py \\
    --reads $reads \\
    --out_dir valid_reads 
    """
}
#!/usr/bin/env nextflow 

process summary_report {
    publishDir "${params.out_dir}", mode: 'copy', pattern: 'report/*.report.tsv'
    label "normal"
    input:
        path(hit_json)
        path(taxon_json)
        path(rpm_json)
    output:
        path("report/*.report.tsv")
    script:
        """
        python $workflow.projectDir/scripts/summary_report.py --hit_json ${hit_json} --taxon_json ${taxon_json} --rpm_json ${rpm_json}  --out_dir report
        """
}
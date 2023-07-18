#!/usr/bin/env nextflow 

process abundance_calculation {
    publishDir "${params.out_dir}/rpm", mode: 'copy', pattern: '*.rpm.json'

    input:
        path(hit_json)
        path(taxon_json)
    output:
        path(hit_json), emit: hit_json
        path(taxon_json), emit: taxon_json
        path("rpm/*.rpm.json"), emit: rpm_json
    script:
        """
        python $workflow.projectDir/scripts/abundance_calculation.py --taxon_json ${taxon_json} --out_dir rpm
        """
}

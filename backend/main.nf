#!/usr/bin/env nextflow 
nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    Usage:
    Pathogen detection pipeline
    Options
    --reads Path to glob pattern of reads
    --out_dir Path to output directory
    --threads Number of threads to use
    """.stripIndent()
    
}

if (params.help) {
    helpMessage()
    System.exit(0)
}
process echo_files {
    publishDir params.out_dir, mode: 'copy', pattern: '*.txt'
    input:
        file(reads)
    script:
        """
        echo $reads > reads.txt
        """
}
include { validate_reads } from "./workflows/modules/validate.nf"
include { host_filter; subsample_reads } from "./workflows/modules/host_filter.nf"
include { pathogen_alignment; assign_taxon } from "./workflows/modules/pathogen_alignment.nf"
include { abundance_calculation } from "./workflows/modules/abundance_calculation.nf"
include { summary_report } from "./workflows/modules/summary_report.nf"  

script_dir = file("$workflow.projectDir").getParent() 
println "Project : $workflow.projectDir"
println "scripts: $workflow.projectDir/scripts"
println "reads: $params.reads"
println "output: $params.out_dir"
println "threads: $params.threads"

workflow {

    reads = channel.fromPath(params.reads, checkIfExists: true)
        .map {it}
        .view { "reads: $it" }

    validated_reads = validate_reads(reads)
    
    nonhuman_ch = host_filter(validated_reads.valid_reads)
    subsampled_ch = subsample_reads(nonhuman_ch.nonhuman_reads)

    hit_ch = pathogen_alignment(subsampled_ch.subsampled_reads.collect()).flatten()
    taxon_ch = assign_taxon(hit_ch)
    abundance_ch = abundance_calculation(taxon_ch.hit_json, taxon_ch.taxon_json)
    summary_report(abundance_ch.hit_json, abundance_ch.taxon_json, abundance_ch.rpm_json)
        .view { "summary report: $it" }

}
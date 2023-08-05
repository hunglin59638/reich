#!/usr/bin/env nextflow 

process host_filter {
    label "performance"
    cpus = 3
    publishDir "${params.out_dir}", mode: 'copy', pattern: 'nonhuman/*.nonhuman.fq'
    input:
        path(reads)    
    output:
        path("nonhuman/*.qc.nonhuman.fq"), emit: nonhuman_reads
    
    script:
    """
    python $workflow.projectDir/scripts/host_filter.py \\
    --fastq_1 $reads \\
    --out_dir nonhuman \\
    --bowtie2_idx $params.bowtie2_idx \\
    --hisat2_idx $params.hisat2_idx \\
    --threads ${task.cpus} 
    """

}

process subsample_reads {
    label "normal"
    publishDir "${params.out_dir}", mode: 'copy', pattern: 'subsampled_reads/*.subsampled.fq'
    input:
        path(reads)
    output:
        path("subsampled_reads/*.subsampled.fq"), emit: subsampled_reads
    script:
    """
    python $workflow.projectDir/scripts/subsample.py \\
    --reads $reads \\
    --out_dir subsampled_reads
    """
} 
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FastQC {
    label 'fastqc'
    
    input:
    path data
    
    output:
    path "results/*" 
    
    script:
    """
    mkdir -p results
    fastqc --noextract -o results ${data}
    """
}

process MultiQC {
    label 'multiqc'
  
    input:
    path results_files
    
    output: 
    path "multiqc_report.html" 

    script:
    """
    multiqc .
    """
}

process trimmomatic {

    label 'trimmomatic'

    input:
    tuple val(name), path(reads)

    output:
    path("${name}_trimmed.fq.gz"), emit: trimmed_reads

    script:
    """
    trimmomatic SE -phred33 $reads ${name}_trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}


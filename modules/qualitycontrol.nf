#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Fastqc_before {
    label 'fastqc'
    
    input:
    tuple val(name), path(reads)
    
    output:
    path "*_fastqc.zip"
    
    script:
    """
    mkdir -p results
    fastqc --noextract -o ${reads}
    """
}

process Multiqc_before {
        
        PublishDir = "$baseDir/results/qc"
        
        label 'multiqc'
    
        input:
        path "*_fastqc.zip"
    
        output: 
        path "raw_multiqc_report.html" 

        script:
        """
        multiqc .
        """
        }
        
        
process Fastqc_after{
        
        label 'fastqc'
        
        input:
        tuple val(name), path(trimmed_reads_1)
        tuple val(name), path(trimmed_reads_2)
        
        output:
        path "*_fastqc.zip"
        
        script:
        """
        fastqc --noextract -o $trimmed_reads_1 $trimmed_reads_2
        """
        }
        
process Multiqc_after{
        
        PublishDir = "$baseDir/results/qc"
        
        label 'multiqc'
        
        input:
        path "*_fastqc.zip"
        
        output:
        path "processed_multiqc.report.html"
        
        script:
        """
        multiqc . "processed_multiqc"
        """
        }
        
        
process trimmomatic {

    publishDir "${baseDir}/results", pattern: "*_trimlog.txt", saveAs: { filename -> filename }

    label 'trimmomatic'

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("*_val_1_fq.gz"), emit: trimmed_reads_1
    tuple val(name), path("*_val_2_fq.gz"), emit: trimmed_reads_2

    script:
    """
    trimmomatic PE -phred33 $reads[0] $reads[1] ${name}_val_1.fq.gz ${name}_unpaired_1.fq.gz ${name}_val_2.fq.gz ${name}_unpaired_2.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
    }


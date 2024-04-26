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

process trim_galore {

    label 'trim_galore'

    input:
    tuple val(name), path(reads)

    output:
     tuple val(name), path("${name}_trimmed.fq.gz"), emit: trimmed_reads

    script:
    """
    trim_galore --quality 20 --phred33 ${reads} --basename ${name} --small_rna 
    """
  
}
//--path_to_cutadapt

process fastqc_trimmed {
    label 'fastqc'
    
    input:
    path trimmed_reads
    
    output:
    path "results/*" 
    
    script:
    """
    mkdir -p results
    fastqc --noextract -o results ${data}
    """
}

process make_transposable_element_gene {
    label 'seqkit'
    
    input:
    path(latest_transcriptome)
    path(Araport11)
    
    output:
    path "transposable_elements_atRTD3.fasta", emit: ref
    
    script:
    """
    grep transposable_element_gene ${Araport11} | cut -f 1 > filter_these_transposable_element_gene.txt
    
    seqkit grep -r -f filter_these_transposable_element_gene.txt ${latest_transcriptome} -o transposable_elements_atRTD3.fasta 
    """
    
}


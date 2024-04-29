#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FastQC {
    label 'fastqc'
    
    input:
      tuple val(name), path(data)
    
    output:
      path("*")
    
    script:
    """
  
    fastqc --noextract -o . ${data}
    """
}

process MultiQC {
    label 'multiqc'
    publishDir "results/multiqc"

    input:
      path(results_files)
    
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
     tuple val(name), path("*trimmed.fq.gz"), emit: trimmed_reads

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

process bowtie_index {
  label 'bowtie'
  
  input:
  path(ref)
  
  output:
  path "*"
  
  script:
  """
  bowtie-build ${ref} ${ref.getSimpleName()}
  """

}

process bowtie_align {
  label 'bowtie'
  
  input:
  tuple val(name), path(trimmed_reads)
  path(index)
  
  output:
  tuple val(name), path("*.fastq"), emit: filtered_reads
  tuple val(name), path("*sam"), emit: sam
  path("*.log"), emit: logs
  
  script:
  """
  bowtie -x ${index[0].getSimpleName()} ${trimmed_reads} -q --un ${trimmed_reads.getSimpleName()}_filtered.fastq -S ${name}.sam 2> ${name}.log
  """
}

process sam_to_bam {
  label 'samtools'
  
  input:
  tuple val(name), path(alignment)
  
  output:
  tuple val(name), path("*.bam"), emit: bam
  
  script:
  """
  samtools view -b -o ${name}.bam ${alignment}
  """
}
  
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
    publishDir "results_TAIR10/trim_galore/", pattern: "*.html"
    
    input:
    tuple val(name), path(reads)

    output:
     tuple val(name), path("*trimmed.fq.gz"), emit: trimmed_reads
     tuple val(name), path("*trimmed_fastqc.html"), emit: trimmed_fastqc
     
    script:
    """
    trim_galore --quality 20 --phred33 ${reads} --basename ${name} --fastqc --illumina
    """
  
}
//--path_to_cutadapt

//process fastqc_after {
    //label 'fastqc'
    
    //input:
     //tuple val(name), path(data)
    
    //output:
    //path "*" 
    
    //script:
    //"""
    //fastqc --noextract -o ${data}
    //"""
//}

//process multiqc_after {
    //label 'multiqc'
    //publishDir "results/multiqc_after"

    //input:
     //path "*_fastqc.zip"     
    
    //output: 
     //path "multiqc_report.html" 

    //script:
    //"""
    //multiqc . -n "after_processing_multiqc_report.html"
    //"""
//}


//process make_transposable_element_gene {
    //label 'seqkit'
    
    //input:
    //path(latest_transcriptome)
    //path(Araport11)
    
    //output:
    //path "transposable_elements_atRTD3.fasta", emit: ref
    
    //script:
    //"""
    //grep transposable_element_gene ${Araport11} | cut -f 1 > filter_these_transposable_element_gene.txt
    
    //seqkit grep -r -f filter_these_transposable_element_gene.txt ${latest_transcriptome} -o transposable_elements_atRTD3.fasta 
    //"""
    
//}


process bowtie_index {
    label 'bowtie'
    publishDir "results_TAIR10/index", pattern: "*.ebwt", mode: 'copy'
  
    input:
     path(ref)
  
    output:
    tuple path("*.ebwt"), val(ref.baseName), emit: index_files
  
    script:
    """
    bowtie-build ${ref} ${ref.baseName}
    """
}

process bowtie_align {
    label 'bowtie'
    publishDir "results_TAIR10/filtering/", pattern: "*.sam", mode: 'copy'
  
    input:
    tuple val(name), path(trimmed_reads)
    tuple path(index_files), val(ref_name)

    output:
    tuple val(name), path("*.fastq"), emit: filtered_reads
    tuple val(name), path("*sam"), emit: sam
    path("*.log"), emit: logs
  
    script:
    """
    bowtie -x ${ref_name} ${trimmed_reads} -q --un ${trimmed_reads.getSimpleName()}_filtered.fastq -S ${name}.sam 2> ${name}.log
    """
}


process sam_to_bam {
  label 'samtools'
  publishDir "results_TAIR10/bam/", pattern: "*.bam", mode: 'copy'
  
  input:
  tuple val(name), path(alignment)
  
  output:
  tuple val(name), path("*.bam"), emit: bam
  
  script:
  """
  samtools view -b -o ${name}.bam ${alignment}
  """
}

//process feature_counts {
  
  //label 'subread'
  //publishDir "results/feature_counts", pattern: "*.txt"

  //input:
  //tuple val(name), path(bam)

  //output: 
  //tuple val(name), path("*.txt"), emit: txt

  //script:
  //"""
  //samtools sort -o ${name}_sorted.bam ${bam}
  //featureCounts -T 5 -t exon -g gene_id -a annotation.gtf -o ${name}_counts.txt ${name}_sorted.bam
  //""" 
//}

//process feature_counts {
    //label 'subread'
    //publishDir "results/feature_counts", pattern: "*.txt"

    //input:
    //tuple val(name), path(bam)

    //output: 
    //tuple val(name), path("*.txt"), emit: txt

    //script:
    //"""
    //featureCounts -T 5 -t exon -g gene_id -a /vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/siRNA_annotate.gff -o ${name}_counts.txt ${bam}
    //"""
//}


process bam_to_sorted_bam {
  
  label 'samtools'
  publishDir "results_TAIR10/sorted_bam/", pattern: "*.sorted.bam", mode: 'copy'
  
  input:
  tuple val(name), path(bam_file)
  
  output:
  tuple val(name), path("*sorted.bam"), emit: sorted_bam
  
  script:
  """
  samtools sort -o ${name}.sorted.bam ${bam_file}
  """
}  

process sorted_bam_to_index {
  
    label 'samtools'
    publishDir "results_TAIR10/sorted_bam/", pattern: "*.sorted.bam.bai", mode: 'copy'
    
    input:
    tuple val(name), path(sorted_bam_file) 
    
    output:
    tuple val(name), path("${sorted_bam_file}.bai"), emit: bam_index
    
    script:
    """
    samtools index ${sorted_bam_file}
    """
}

process idxstats {
    
    label 'samtools'
    publishDir "results_TAIR10/idxstats/", pattern: "*.log", mode: 'copy'
    
    input:
    tuple val(name), path(sorted_bam_file)
    
    output:
    tuple val(name), path("${sorted_bam_file}.log"), emit: bam_log
    
    script:
    """
    samtools idxstats ${sorted_bam_file} > ${sorted_bam_file}.log
    """
}

//process sorted_bam_to_bed {

   //label 'bedtools'
   
   //input:
   //tuple val(name), path(sorted_bam_file)

   //output:
   //tuple val(name), path("${sorted_bam_file}.bed"), emit: bam.bed

   //script:
   //"""
   //bedtools bamtobed -i ${sorted_bam_file} > ${sorted_bam_fiel}.bed
   //"""
//}
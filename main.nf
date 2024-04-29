#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FastQC}                         from "./modules/qualitycontrol.nf"
include {MultiQC}                        from "./modules/qualitycontrol.nf"
include {trim_galore}                    from "./modules/qualitycontrol.nf"
include {fastqc_trimmed}                 from "./modules/qualitycontrol.nf"
include {make_transposable_element_gene} from "./modules/qualitycontrol.nf"
include {bowtie_index}                   from "./modules/qualitycontrol.nf"
include {bowtie_align}                   from "./modules/qualitycontrol.nf"
include {sam_to_bam}                     from "./modules/qualitycontrol.nf"

workflow {
    
    
    
    myFileChannel = Channel.fromFilePairs('raw_data/*/*/*.fq.gz', size: 1, flat: true)
    myFileChannel.view()
    FastQC(myFileChannel)
    MultiQC(FastQC.out.collect())
    trim_galore(myFileChannel)
    //multiqc_after(trim_galore.out)
    
    
    
    latest_transcriptome = Channel.fromPath("trancriptome.fa")
    araport11 = Channel.fromPath("Araport11_gene_type")
    make_transposable_element_gene(latest_transcriptome, araport11)

    bowtie_index(make_transposable_element_gene.out.ref)
    bowtie_align(trim_galore.out.trimmed_reads, bowtie_index.out.collect())
    
    sam_to_bam(bowtie_align.out.sam)
    
}    
    //FastQC_results = FastQC(myFileChannel).collect()
    //MultiQC(FastQC_results)  // Use the specific output channel from FastQC
    //trim_galoreChannel = Channel.fromPath('raw_data/*/*/*.fq.gz')
    //trim_galoreChannel = Channel.fromFilePairs('raw_data/*/*/*.fq.gz', size: 1, flat: true)
    //trim_galoreChannel.view()
    //trimmed_reads = trim_galore(trim_galoreChannel)
    
    
    //make_transposable_element_gene(latest_transcriptome, araport11)
    //bowtie_index(make_transposable_element_gene.out.ref)
    //bowtie_align(trimmed_reads.out, bowtie_index.out.collect())
    
    //FastQC_results_trimmed = FastQC(trimmed_reads)


//trimmomatic.out.trimmed_reads


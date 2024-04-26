#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FastQC} from "./modules/qualitycontrol.nf"
include {MultiQC} from "./modules/qualitycontrol.nf"
include {trim_galore} from "./modules/qualitycontrol.nf"
include {fastqc_trimmed} from "./modules/qualitycontrol.nf"
include {make_transposable_element_gene} from "./modules/qualitycontrol.nf"
include {bowtie_index} from "./modules/qualitycontrol.nf"

workflow {
  
    
    //myFileChannel = Channel.fromPath('raw_data/*/*/*_[1,2,3].fq.gz')
    //FastQC_results = FastQC(myFileChannel).collect()
    //MultiQC(FastQC_results)  // Use the specific output channel from FastQC
    //trim_galoreChannel = Channel.fromPath('raw_data/*/*/*.fq.gz')
    //trim_galoreChannel = Channel.fromFilePairs('raw_data/*/*/*.fq.gz', size: 1, flat: true)
    //trim_galoreChannel.view()
    //trimmed_reads = trim_galore(trim_galoreChannel)
    
    latest_transcriptome = Channel.fromPath("trancriptome.fa")
    araport11 = Channel.fromPath("Araport11_gene_type")
    make_transposable_element_gene(latest_transcriptome, araport11)
    
    bowtie_index(make_transposable_element_gene.out.ref)
    
    
    //FastQC_results_trimmed = FastQC(trimmed_reads)
}

//trimmomatic.out.trimmed_reads


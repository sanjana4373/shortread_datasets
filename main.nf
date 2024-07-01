#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FastQC}                         from "./modules/qualitycontrol.nf"
include {MultiQC}                        from "./modules/qualitycontrol.nf"
include {trim_galore}                    from "./modules/qualitycontrol.nf"
//include {fastqc_after}                   from "./modules/qualitycontrol.nf"
//include {multiqc_after}                  from "./modules/qualitycontrol.nf"
//include {make_transposable_element_gene} from "./modules/qualitycontrol.nf"
include {bowtie_index}                   from "./modules/qualitycontrol.nf"
include {bowtie_align}                   from "./modules/qualitycontrol.nf"
include {sam_to_bam}                     from "./modules/qualitycontrol.nf"
include {bam_to_sorted_bam}              from "./modules/qualitycontrol.nf"
include {sorted_bam_to_index}            from "./modules/qualitycontrol.nf"
include {idxstats}                       from "./modules/qualitycontrol.nf"

//include {feature_counts}                 from "./modules/qualitycontrol.nf"
//include {bowtie2_index}                    from "./modules/qualitycontrol.nf"
//include {bowtie2_align}                    from "./modules/qualitycontrol.nf"
//include {sorted_bam_to_bed}              from "./modules/qualitycontrol.nf"


workflow {
    
    myFileChannel = Channel.fromFilePairs('raw_data/*/*/*.fq.gz', size: 1, flat: true)
    myFileChannel.view()
    FastQC(myFileChannel)
    MultiQC(FastQC.out.collect())
    trim_galore(myFileChannel)
    //fastqc_after(trim_galore.out)
    //multiqc_after(trim_galore.out.trimmed_fastqc.collect())
    
    whole_genome = Channel.fromPath('/vol/*/*/*all.fasta')
    //whole_genome = Channel.fromPath('/vol/*/*/*/*/reference_genome/*.fa') 
    //latest_transcriptome = Channel.fromPath("trancriptome.fa")
    //araport11 = Channel.fromPath("Araport11_gene_type")
    //make_transposable_element_gene(latest_transcriptome, araport11)
    
    bowtie_index(whole_genome)
    bowtie_align(trim_galore.out.trimmed_reads, bowtie_index.out.collect())
   
    sam_to_bam(bowtie_align.out.sam)
    bam_to_sorted_bam(sam_to_bam.out.bam)

    sorted_bam_to_index(bam_to_sorted_bam.out.sorted_bam)
    idxstats(bam_to_sorted_bam.out)
    
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


    sorted_bam_to_index(bam_to_sorted_bam.out.sorted_bam)
    idxstats(bam_to_sorted_bam.out)
    //feature_counts(sam_to_bam.out.bam)
    //feature_counts(bam_to_sorted_bam.out.sorted_bam)

} 


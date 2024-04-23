#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FastQC} from "./modules/qualitycontrol.nf"
include {MultiQC} from "./modules/qualitycontrol.nf"
include {trimmomatic} from "./modules/qualitycontrol.nf"


workflow {
    trimmomaticChannel = Channel.fromFilePairs('raw_data/*/*/*.fq.gz', size: 1, flat: true)
    myFileChannel = Channel.fromPath('raw_data/*/*/*.fq.gz')
    FastQC_results = FastQC(myFileChannel).collect()
    MultiQC(FastQC_results)  // Use the specific output channel from FastQC
    trimmomatic(trimmomaticChannel)
}
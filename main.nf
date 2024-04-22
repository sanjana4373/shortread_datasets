main.nf
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

workflow {
    myFileChannel = Channel.fromPath('raw_data/*/*/*.fq.gz')
    FastQC_results = FastQC(myFileChannel).collect()
    MultiQC(FastQC_results)  // Use the specific output channel from FastQC
}


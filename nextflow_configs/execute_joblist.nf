#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.joblist = 'NONE'

process execute_jobs {

    errorStrategy 'retry'
    maxRetries {}

    input:
    val line

    // one line represents an independent command
    script:
    """
    ${line}
    """
}

workflow {
    if (params.joblist == "NONE"){
        println("Usage: nextflow execute_joblist.nf  --joblist [joblist file] -c [config file]")
        System.exit(2);
    }

    def lines = Channel.fromPath(params.joblist).splitText()
    execute_jobs(lines)
}

#!/usr/bin/env nextflow

samples = '*.bam*'

process svaba_calling {
	cpus 2	
	maxForks 8
	module 'svaba'
	publishDir "svaba_vcf"
	input:
		tuple val(core), path(f)
	output:
		path("${core}.vcf.gz")
	script:
	out="${core}.vcf.gz"
	"""
	mkdir ${core}.vcf.gz
	chmod a+wrx ${core}.vcf.gz
	svaba run -p 4 -G ${params.genome} -t ${f[0]} -z ./
	"""	
 }

workflow {
    input = Channel.fromFilePairs(params.input+"/"+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		svaba_calling(input)
	        input.view()

}
                                              

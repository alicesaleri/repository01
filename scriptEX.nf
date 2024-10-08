#!/usr/bin/env nextflow

samples = '*.bam*'

process insurveyor_calling {
	cpus 4	
	maxForks 8
	module 'INSurVeyor'
	input:
		tuple val(core), path(f)
	output:
		path("${core}.vcf.gz")
	script:
	out="${core}.vcf.gz"
	"""
	insurveyor.py -t 4 ${f[0]} ./ ${params.genome}
	mv out.vcf.gz ${out}
	"""
 }

workflow {
    input = Channel.fromFilePairs(params.input+"/"+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		insurveyor_calling(input)
	        input.view()


}
                                              

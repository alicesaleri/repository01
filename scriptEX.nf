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
	echo "running!!!!!!"
	insurveyor.py --threads 4 ${f[0]} ./ ${params.genome}
	echo "insurveyor.py completed!!!!!!!!!"
	"""
 }

workflow {
    input = Channel.fromFilePairs(params.input+"/"+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		insurveyor_calling(input)
	        /*input.view()*/


}
                                              

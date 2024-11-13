#!/usr/bin/env nextflow

samples = '*.bam*'
ensamples = '*.vcf*'
//samples = 'LP6005443-DNA_D01*bam*'

process insurveyor_calling {
	cpus 8
	maxForks 4
	module 'INSurVeyor'
	memory '40GB'
	publishDir "insurveyor_vcf"
	input:
		tuple val(core), path(f)
	output:
		path("${core}.insurveyor.vcf.gz")
		path("${core}.insurveyor.vcf.gz.tbi")
	script:
	out="${core}.insurveyor.vcf.gz"
	"""
	set -euxo pipefail
	insurveyor.py --threads 8 ${f[0]} ./ ${params.genome}
	mv out.pass.vcf.gz ${out}
	tabix -p vcf ${out}
	"""	
 }

process svaba_calling {
	errorStrategy 'finish'
	cpus 24
	maxForks 4
	module 'svaba'
	memory '40GB'
	publishDir "svaba_vcf"
	input:
		tuple val(core), path(f)
	output:
  		path("${core}.svaba.vcf.gz")
		path("${core}.svaba.vcf.gz.tbi")  
	script:
	out1="${core}.svaba.sv.vcf.gz"
	out2="${core}.svaba.indel.vcf.gz"
	out3="${core}.svaba.vcf.gz"
	"""
	set -euxo pipefail
	/usr/bin/time -o "resources.t" -f "%e %M" svaba run -p 24 -G ${params.genome} -t ${f[0]} -z 
	mv no_id.svaba.sv.vcf.gz ${out1}
	mv no_id.svaba.indel.vcf.gz ${out2}
	bcftools index -t ${out1}
	bcftools index -t ${out2}
	bcftools concat ${out1} ${out2} -a -D -Oz -o ${out3}
	tabix -p vcf ${out3} 
	"""	
 }

process fix_files{
	publishDir "svaba_fixed_vcf"
	input:
		tuple val (core), path(fs)
	output:
  		path(result)
		path("${result}.tbi")  
	script:
	result="${core}.svaba.fixed.vcf.gz"
	
	"""
	set -euxo pipefail
	zcat "${core}.svaba.vcf.gz" | sd "ID=GQ,Number=1,Type=String" "ID=GQ,Number=1,Type=Integer" | sed 's/ID=PL,Number=[^,]*/ID=PL,Number=G/' | bgzip > ${result}
	tabix -p vcf ${result}	
	"""
}

//	zcat "${core}.svaba.vcf.gz" | sd "ID=GQ,Number=1,Type=String" "ID=GQ,Number=1,Type=Float" | bgzip > ${result}


process truvari_ensemble {
	cpus 2
	maxForks 4
	module 'truvari'
	memory '20GB'
	publishDir "truvari_vcf"
	input:
                path(files)
	output:
		path("tru_coll_sorted.vcf.gz")

		path("tru_coll_sorted.vcf.gz.tbi")
	script:
	"""
	set -euxo pipefail
	
	vcf_files=""
	for file in ${files}/*.vcf.gz; do
		vcf_files="\$vcf_files \$file"
	done
	bcftools merge -m none \$vcf_files -Oz -o bcftools.merge.vcf.gz --force-samples
	bcftools index -t bcftools.merge.vcf.gz
 	truvari collapse -i bcftools.merge.vcf.gz -o truvari.vcf.gz -c tru_collapse.vcf.gz
	vcf-sort tru_collapse.vcf.gz | bgzip -c > tru_coll_sorted.vcf.gz
	tabix -p vcf tru_coll_sorted.vcf.gz 
	"""	
 }
//truvari bench --reference ${params.genome} --base ${files[0]} --comp ${files[2]} --comp ${files[4]} --comp ${files[6]} --comp ${files[8]} -o truBench 
	
//	tabix -p vcf ${core}.merge.vcf.gz 
//bcftools index -t ${files}/*.vcf.gz

process survivor_ensemble {
	cpus 2
	maxForks 4
	memory '20GB'
	publishDir "survivor_vcf"
	input:
                tuple val(core), path(files)
	output:
		path("${core}.surv_rem.vcf.gz")
		path("${core}.surv_rem.vcf.gz.tbi")
	script:
	"""
	set -euxo pipefail
	zcat "${files[0]}" > vcf01.vcf
	zcat "${files[2]}" > vcf02.vcf
	zcat "${files[4]}" > vcf03.vcf
	zcat "${files[6]}" > vcf04.vcf
	zcat "${files[8]}" > vcf05.vcf

	echo "vcf01.vcf" > file.txt
	echo "vcf02.vcf" >> file.txt	
	echo "vcf03.vcf" >> file.txt
	echo "vcf04.vcf" >> file.txt
	echo "vcf05.vcf" >> file.txt

	
	SURVIVOR merge file.txt 1000 1 1 -1 -1 -1 ${core}.survivor.vcf
	bgzip ${core}.survivor.vcf
	zgrep -v HLA ${core}.survivor.vcf.gz | bgzip |  bcftools sort  -Oz -o ${core}.survivor_sorted.vcf.gz
	bcftools view -e 'INFO/END < POS' ${core}.survivor_sorted.vcf.gz -o ${core}.surv_rem.vcf.gz
	bcftools index -t ${core}.surv_rem.vcf.gz
	"""	
}
//bcftools norm -c ws -f ${params.genome} ${core}.survivor_sorted.vcf.gz -o ${core}.surv_norm.vcf.gz

process survclusterer_ensemble {
	cpus 2
	maxForks 4
	memory '20GB'
	publishDir "clusterer_vcf"
	input:
                tuple val(core), path(files)
	output:
		path("${core}.clusterer.vcf.sv.gz")
		path("${core}.clusterer.vcf.sv.gz.tbi")
	script:
	"""
	set -euxo pipefail
	echo "${files[0]}" > file.txt
	echo "${files[2]}" >> file.txt	
	echo "${files[4]}" >> file.txt
	echo "${files[6]}" >> file.txt
	echo "${files[8]}" >> file.txt

	clusterer -d 1000 -t 4 file.txt ${params.genome} -o ${core}.clusterer.vcf
	bgzip -k ${core}.clusterer.vcf.sv
	tabix ${core}.clusterer.vcf.sv.gz
	"""	
 }

process annotsv_annotation {
	cpus 2
	maxForks 4
	module 'AnnotSV/3.4'
	memory '40GB'
	publishDir "annotsv_vcf"
	input:
                path(files)
	output:
		path("AnnotSV_dir/*")
	script:
	"""
	set -euxo pipefail

	gunzip -k -c ${files}/*coll*.vcf.gz > tru.vcf 
	mkdir -p AnnotSV_dir
	${ANNOTSV_DIR} -genomeBuild  GRCh38 -SVinputFile tru.vcf -outputDir AnnotSV_dir -outputFile annotsv -SVinputInfo 1
	"""	
 }


workflow {
	input = Channel.fromFilePairs(params.input+samples)
   		{ fname -> fname.simpleName.replaceAll(".md.*", "")}
	input2 = Channel.fromPath("/home/saleri/project_cnv_vr/scriptNF_running/svaba_vcf/*.vcf.gz").map {f -> 
		def core= f.getBaseName().toString().split("\\.")[0] 
		return tuple(core,f)}
	ensinput = Channel.fromFilePairs(params.directories.collect{it+ensamples}, size: 10)
{fname -> fname.simpleName.replaceAll("","")}
	merinput = Channel.fromPath("/home/saleri/project_cnv_vr/scriptNF_running/survivor_vcf/")	
	anninput = Channel.fromPath("/home/saleri/project_cnv_vr/scriptNF_running/truvari_vcf")

	main:
		//insurveyor_calling(input)
		//svaba_calling(input)
		//fix_files(input2)
		//survivor_ensemble(ensinput)
		//survclusterer_ensemble(ensinput)
		//truvari_ensemble(merinput)
		annotsv_annotation(anninput)		
		//println(params.directories.collect{it+ensamples})
		//anninput.view()


}
                                       
                                              
//commit 13th november 

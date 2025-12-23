#!/usr/bin/env nextflow

/* Organize summary stats for fGWAS */

nextflow.enable.dsl=2
def channel_glob(p) {
	chan = Channel.fromPath(p[0], type: 'any').map{it -> [it.name.replaceAll(p[1], ""), it]}
	return chan
}

def channel_glob_split(p, delim) {
	chan = Channel.fromPath(p[0], type: 'any')
	    .map{it -> it.name.replaceAll(p[1], "").tokenize(delim).plus([it])}
	return chan
}

def getfile(p){
	return Channel.fromPath(p)
}

def abspath(f){
	File file = new File(f)
	path = file.getAbsolutePath()
	return path
}

workflow {
	gwas = channel_glob(params.gwas)
	munge_stats(gwas) // trait summary

	annots = channel_glob(params.annotations) // annot annotbed
	intersect_annots(munge_stats.out.tointersect.combine(annots)) // trait traitbed annot annotbed

	collapse_in = intersect_annots.out.groupTuple(by: 0).combine(munge_stats.out.main, by: 0)
	collapse_intersects(collapse_in)	

	annotlist = annots.map{it -> [it[0]]}
	fgwas_single_annot(collapse_intersects.out.combine(annotlist))
	plot(fgwas_single_annot.out.main.groupTuple(by: 0))
}

process munge_stats {

	publishDir "${params.results}/sumstats"
	memory '20GB'
	errorStrategy "finish"
	time '2h'
	tag "${trait}"

	input:
	tuple val(trait), path(sumstats)

	output:
	tuple val(trait), path("${trait}.sumstats.txt"), emit: main
	tuple val(trait), path("${trait}.sumstats.bed"), emit: tointersect

	"""
	munge_sumstats.py --sumstats $sumstats --prefix ${trait}.sumstats  --N n_complete_samples \
	--frq EAF --effect beta --stderr se --pos snp_end 
	"""

}

process intersect_annots {
	/*Get 0/1 columns for annotation overlap for SNPs */

	storeDir "${params.results}/intersect_temp"
	memory '4GB'
	errorStrategy "finish"
	time '2h'
	tag "${trait}__annot__${annot}"
	
	input:
	tuple val(trait), path(summarybed), val(annot), path(annotbed)

	output:
	tuple val(trait), path("${prefix}.bed")

	script:
	prefix = "${trait}__annot__${annot}"

	"""
	intersectBed -a $summarybed -b $annotbed -c | awk '{ if ((\$NF>1)) {\$NF=1}; print \$NF }' > temp;
	echo -e "$annot" | cat - temp > ${prefix}.bed ;
	rm temp
	"""

}


process collapse_intersects {
	/*Assemble all input into one file */

	publishDir "${params.results}/formatted_input", mode: "rellink"
	memory '4GB'
	errorStrategy "finish"
	time '2h'
	tag "${trait}"

	input:
	tuple val(trait), path(annotbed), path(summary)

	output:
	tuple val(trait), path("${trait}.input.gz")
	
	"""
	paste -d' ' $summary ${annotbed.join(' ')} | gzip -c > ${trait}.input.gz ;
	"""

}


process fgwas_single_annot {
	/*Run fGWAS for each annotation separately */
	
	publishDir "${params.results}/single_annot_fgwas", mode: 'rellink'
	memory '10GB'
	errorStrategy "finish"
	time '10h'
	tag "${trait}__annot__${annot}"
	
	input:
	tuple val(trait), path(input), val(annot) 
	
	output:
	tuple val(trait), path("${prefix}.params"), path("${prefix}.llk"), emit: main
	path("${prefix}.log")

	script:
	prefix = "${trait}__annot__${annot}"
	runtype = (params.type == "eqtl") ? " -fine " : ""
	
	"""
	fgwas  -i $input  -w $annot	 -o ${prefix} ${runtype} &> ${prefix}.log 
	"""

}

process plot {
	/*Compile fGWAS results for selected annotations per iteration */
	publishDir "${params.results}/figures", mode: 'rellink'
	memory '4GB'
	errorStrategy "finish"
	time '2h'
	tag "${trait}"

	input:
	tuple val(trait), path(pars), path(llk)

	output:
	tuple val(trait), path("${trait}.results*") 
	
	"""
	compile_single.py --pars ${pars.join(' ')} --llk ${llk.join(' ')} --trait $trait --output ${trait}.results.tsv ;
	plot_single.R  ${trait}.results.tsv	 ${trait}.results
	"""
	
}

workflow.onComplete {
	if (workflow.success){
		subject = "fGWAS organize execution complete"
	}
	else {
		subject = "fGWAS execution error"
	}

	recipient = params.email

	['mail', '-s', subject, recipient].execute() << """

	Pipeline execution summary
	---------------------------
	Completed at: ${workflow.complete}
	Duration	: ${workflow.duration}
	Success		: ${workflow.success}
	workDir		: ${workflow.workDir}
	exit status : ${workflow.exitStatus}
	Error report: ${workflow.errorReport ?: '-'}
	"""
}

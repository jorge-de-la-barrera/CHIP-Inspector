## Jorge de la Barrera, 2020
##
## Call somatic variants for a single sample

workflow somatic_calling {
	String sample_id
	String queue
	String accounting
	Object wf
	Object genome
	String references_dir
	String vcf_in

	call FilterMutectCalls{
		input:
 			sample_id =sample_id,
			execution = wf.FilterMutectCalls,
			queue = queue,
			accounting = accounting,
			genome = genome,
			references_dir = references_dir,
			vcf_in = vcf_in
	}
	output {
		File filtered_vcf_out = FilterMutectCalls.filtered_vcf_out
		File filtered_vcf_out_idx = FilterMutectCalls.filtered_vcf_out_idx
		File filtering_stats = FilterMutectCalls.filtering_stats
	}
}

task FilterMutectCalls {
	String sample_id
	Object execution
	String queue
	String accounting
	Object genome
	String references_dir
	String vcf_in

	String genome_fasta = references_dir + "/" + genome.path + "/" + genome.ref_fasta

	command <<<
		# Load enviroment
		source /programs/GATK/env_vep.sh

		set -e

		gatk --java-options "${execution.java_args_FilterMutectCalls}" FilterMutectCalls \
			${"-R " + genome_fasta} \
			${"-V " + vcf_in} \
			--filtering-stats "${sample_id}.vcf.gz.filteringStats.tsv" \
			-O "${sample_id}.vcf.gz"

	>>>
	runtime {
		backend : 'SGE_nope_long'
		cpu : execution.cpu
		memory : execution.memory
		accounting : accounting
		sge_project : queue
	}
	output {
		File filtered_vcf_out = "${sample_id}.vcf.gz"
		File filtered_vcf_out_idx = "${sample_id}.vcf.gz.tbi"
		File filtering_stats = "${sample_id}.vcf.gz.filteringStats.tsv"
	}
}


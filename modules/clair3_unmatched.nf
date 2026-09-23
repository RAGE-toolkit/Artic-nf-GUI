// modules/clair3_unmatched.nf
process CLAIR3_UNMATCHED {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_bam), path(reference), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.unmatched.vcf"), emit: vcf

	script:
		if( !input_bam.exists() )
			exit 1, "CLAIR3_UNMATCHED: Input BAM not found: ${input_bam}"
		if( !reference.exists() )
			exit 1, "CLAIR3_UNMATCHED: Reference not found: ${reference}"

		"""
		set -euo pipefail
		[ -s "${input_bam}.bai" ] || samtools index "${input_bam}"

		run_clair3.sh \
			--enable_long_indel \
			--chunk_size=10000 \
			--haploid_sensitive \
			--no_phasing_for_fa \
			--bam_fn="${input_bam}" \
			--ref_fn="${reference}" \
			--output="clair3_unmatched" \
			--threads="${params.threads}" \
			--platform=ont \
			--model_path="${params.model_path}" \
			--include_all_ctgs \
			--enable_variant_calling_at_sequence_head_and_tail

		bgzip -dc clair3_unmatched/merge_output.vcf.gz > "${params.run_name}_${sampleId}.unmatched.vcf"
		"""
}

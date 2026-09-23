// modules/clair3-2.nf
process CLAIR3_2 {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_bam), path(reference), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.2.vcf"), emit: vcf

	script:
		if( !input_bam.exists() )
			exit 1, "CLAIR3_2: Input BAM not found: ${input_bam}"
		if( !reference.exists() )
			exit 1, "CLAIR3_2: Reference not found: ${reference}"

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
			--output="clair3_pool2" \
			--threads="${params.threads}" \
			--platform=ont \
			--model_path="${params.model_path}" \
			--include_all_ctgs \
			--enable_variant_calling_at_sequence_head_and_tail

		bgzip -dc clair3_pool2/merge_output.vcf.gz > "${params.run_name}_${sampleId}.2.vcf"
		"""
}

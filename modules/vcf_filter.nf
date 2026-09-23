// modules/vcf_filter.nf
process VCF_FILTER {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_vcf), path(vcf_filter_script), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.pass.vcf"),   emit: pass_vcf
		tuple val(sampleId), path("${params.run_name}_${sampleId}.fail.vcf"),   emit: fail_vcf
		tuple val(sampleId), path("${params.run_name}_${sampleId}.ignore.vcf"), emit: ignore_vcf

	script:
		if( !input_vcf.exists() )
			exit 1, "VCF_FILTER: Input VCF not found: ${input_vcf}"
		if( !vcf_filter_script.exists() )
			exit 1, "VCF_FILTER: vcf_filter.py not found: ${vcf_filter_script}"

		"""
		set -euo pipefail

		python "${vcf_filter_script}" \
			--min-depth ${params.mask_depth} \
			--min-variant-quality ${params.min_variant_quality} \
			--min-allele-frequency ${params.min_allele_frequency} \
			--min-mask-allele-frequency ${params.min_mask_allele_frequency} \
			--min-frameshift-quality ${params.min_frameshift_quality} \
			--min-minor-allele-count ${params.min_minor_allele_count} \
			"${input_vcf}" \
			"${params.run_name}_${sampleId}.pass.vcf" \
			"${params.run_name}_${sampleId}.fail.vcf" \
			"${params.run_name}_${sampleId}.ignore.vcf"
    """
}

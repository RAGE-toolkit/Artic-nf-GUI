// modules/vcf_filter.nf
process VCF_FILTER {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_vcf), path(vcf_filter_script), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.pass.vcf"), emit: pass_vcf
		tuple val(sampleId), path("${params.run_name}_${sampleId}.fail.vcf"), emit: fail_vcf
		tuple val(sampleId), path("${params.run_name}_${sampleId}.pass.vcf.gz"), emit: pass_vcf_gz
	// errorStrategy 'terminate'   // prefer to fail loudly; switch to 'ignore' if you must

	script:
		//def vcf_filter = file("${projectDir}/scripts/vcf_filter.py")

		if( !input_vcf.exists() )
			exit 1, "VCF_FILTER: Input VCF not found: ${input_vcf}"
		if( !vcf_filter_script.exists() )
			exit 1, "VCF_FILTER: vcf_filter.py not found: ${vcf_filter_script}"

		"""
		set -euo pipefail

		python "${vcf_filter_script}" --medaka \
			"${input_vcf}" \
			"${params.run_name}_${sampleId}.pass.vcf" \
			"${params.run_name}_${sampleId}.fail.vcf"

		bgzip -fc "${params.run_name}_${sampleId}.pass.vcf" > "${params.run_name}_${sampleId}.pass.vcf.gz"
		tabix -p vcf "${params.run_name}_${sampleId}.pass.vcf.gz"
    """
}


// modules/compress_and_index_vcf.nf
process COMPRESS_AND_INDEX_VCF {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(pass_vcf), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.pass.vcf.gz"), emit: vcf_gz
		path "${params.run_name}_${sampleId}.pass.vcf.gz.tbi",                   emit: vcf_gz_tbi

	script:
		if( !pass_vcf.exists() )
			exit 1, "COMPRESS_AND_INDEX_VCF: Input VCF not found: ${pass_vcf}"

		"""
		set -euo pipefail

		bgzip -kf "${pass_vcf}"
		tabix -f -p vcf "${params.run_name}_${sampleId}.pass.vcf.gz"
		"""
}

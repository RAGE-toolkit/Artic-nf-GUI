// modules/bcftools_norm.nf
process BCFTOOLS_NORM {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(preconsensus_fa), path(pass_vcf_gz), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.normalised.vcf.gz"), emit: normalised_vcf
		tuple val(sampleId), path("${params.run_name}_${sampleId}.normalised.vcf.gz.tbi"), emit: normalised_tbi

	script:
		if( !preconsensus_fa.exists() )
			exit 1, "BCFTOOLS_NORM: Preconsensus FASTA not found: ${preconsensus_fa}"
		if( !pass_vcf_gz.exists() )
			exit 1, "BCFTOOLS_NORM: PASS VCF (bgzipped) not found: ${pass_vcf_gz}"

		"""
		set -euo pipefail

		if [[ ! -f "${preconsensus_fa}.fai" ]]; then
			samtools faidx "${preconsensus_fa}"
		fi

		bcftools norm --check-ref x \
			--fasta-ref "${preconsensus_fa}" \
			-O z -o "${params.run_name}_${sampleId}.normalised.vcf.gz" \
			"${pass_vcf_gz}"

		tabix -f -p vcf "${params.run_name}_${sampleId}.normalised.vcf.gz"
		"""
}

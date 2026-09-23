// modules/bcftools_consensus.nf
process BCFTOOLS_CONSENSUS {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(preconsensus_fa), path(normalised_vcf_gz), path(normalised_vcf_tbi), path(coverage_mask), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path ("${params.run_name}_${sampleId}.consensus.fasta"), emit: consensus_fa

	script:
		if( !preconsensus_fa.exists() )   exit 1, "BCFTOOLS_CONSENSUS: Preconsensus FASTA not found: ${preconsensus_fa}"
		if( !normalised_vcf_gz.exists() ) exit 1, "BCFTOOLS_CONSENSUS: Normalised VCF not found: ${normalised_vcf_gz}"
		if( !coverage_mask.exists() )     exit 1, "BCFTOOLS_CONSENSUS: Coverage mask not found: ${coverage_mask}"

		"""
		set -euo pipefail

		if [[ ! -f "${preconsensus_fa}.fai" ]]; then
			samtools faidx "${preconsensus_fa}"
		fi

		# The .tbi is staged alongside the vcf.gz via the input tuple, but guard
		# against it going missing (e.g. a stale -resume cache entry) since
		# bcftools consensus - unlike samtools - won't regenerate it itself.
		if [[ ! -f "${normalised_vcf_gz}.tbi" ]]; then
			tabix -f -p vcf "${normalised_vcf_gz}"
		fi

		bcftools consensus \
			-f "${preconsensus_fa}" \
			-m "${coverage_mask}" \
			"${normalised_vcf_gz}" \
			-o "${params.run_name}_${sampleId}.consensus.fasta"
		"""
}

// modules/split_unmatched.nf
process SPLIT_UNMATCHED {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(primertrimmed_bam), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.unmatched.primertrimmed.rg.sorted.bam"), emit: unmatched_bam
		path "${params.run_name}_${sampleId}.unmatched.primertrimmed.rg.sorted.bam.bai",                   emit: unmatched_bai

	script:
		if( !primertrimmed_bam.exists() )
			exit 1, "SPLIT_UNMATCHED: Input BAM not found: ${primertrimmed_bam}"

		"""
		set -euo pipefail

		samtools view -b -r unmatched "${primertrimmed_bam}" \
			-o "${params.run_name}_${sampleId}.unmatched.primertrimmed.rg.sorted.bam"
		samtools index "${params.run_name}_${sampleId}.unmatched.primertrimmed.rg.sorted.bam"
		"""
}

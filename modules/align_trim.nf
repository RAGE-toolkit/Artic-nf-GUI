// modules/align_trim.nf
process ALIGN_TRIM {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_bam), val(item), val(scheme), val(version), path(bed)

	output:
		path "${params.run_name}_${sampleId}.alignreport.txt",                                   emit: align_report
		tuple val(sampleId), path("${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam"), emit: primertrimmed_bam
		path "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam.bai",                   emit: primertrimmed_bai
		tuple val(sampleId), path("${params.run_name}_${sampleId}.1.primertrimmed.rg.sorted.bam"), emit: pool1_bam
		path "${params.run_name}_${sampleId}.1.primertrimmed.rg.sorted.bam.bai",                 emit: pool1_bai
		tuple val(sampleId), path("${params.run_name}_${sampleId}.2.primertrimmed.rg.sorted.bam"), emit: pool2_bam
		path "${params.run_name}_${sampleId}.2.primertrimmed.rg.sorted.bam.bai",                 emit: pool2_bai

	script:
		if( !bed.exists() )
			exit 1, "ALIGN_TRIM: Scheme BED not found: ${bed}"
		if( !input_bam.exists() )
			exit 1, "ALIGN_TRIM: Input BAM not found: ${input_bam}"

		"""
		set -euo pipefail

		align_trim --normalise ${params.normalise} "${bed}" \
			--primer-match-threshold ${params.primer_match_threshold} \
			--min-mapq ${params.min_mapq} \
			--allow-incorrect-pairs \
			--report "${params.run_name}_${sampleId}.alignreport.txt" \
			--amp-depth-report "${params.run_name}_${sampleId}.amplicon_depths.tsv" \
			--genome-coverage-report "${params.run_name}_${sampleId}" \
			--samfile "${input_bam}" \
			-o "${params.run_name}_${sampleId}.primertrimmed.rg.bam"

		samtools sort -T "${sampleId}" "${params.run_name}_${sampleId}.primertrimmed.rg.bam" \
			-o "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam"
		rm "${params.run_name}_${sampleId}.primertrimmed.rg.bam"
		samtools index "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam"

		samtools view -b -r 1 "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam" \
			-o "${params.run_name}_${sampleId}.1.primertrimmed.rg.sorted.bam"
		samtools index "${params.run_name}_${sampleId}.1.primertrimmed.rg.sorted.bam"

		samtools view -b -r 2 "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam" \
			-o "${params.run_name}_${sampleId}.2.primertrimmed.rg.sorted.bam"
		samtools index "${params.run_name}_${sampleId}.2.primertrimmed.rg.sorted.bam"
		"""
}

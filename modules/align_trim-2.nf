// modules/align_trim-2.nf
process ALIGN_TRIM_2 {
	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_bam), val(item), val(scheme), val(version), path(bed), path(align_trim)

	output:
		path "${params.run_name}_${sampleId}.alignreport_2.txt",                                      emit: align_report
		tuple val(sampleId), path("${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam"),     emit: primertrimmed_bam
		path "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam.bai",                       emit: primertrimmed_bai

	errorStrategy 'terminate'

	script:
		if( !bed.exists() )
			exit 1, "ALIGN_TRIM_2: Scheme BED not found: ${bed}"
		if( !align_trim.exists() )
			exit 1, "ALIGN_TRIM_2: align_trim.py not found: ${align_trim}"
		if( !input_bam.exists() )
			exit 1, "ALIGN_TRIM_2: Input BAM not found: ${input_bam}"

		"""
		set -euo pipefail

		python "${align_trim}" \
			--normalise ${params.medaka_normalise} \
			"${bed}" \
			--report "${params.run_name}_${sampleId}.alignreport_2.txt" \
			< "${input_bam}" \
			2> "${params.run_name}_${sampleId}.alignreport.err" \
			| samtools sort -T "${sampleId}" -o "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam"

		samtools index "${params.run_name}_${sampleId}.primertrimmed.rg.sorted.bam"
    """
}


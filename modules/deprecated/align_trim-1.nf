// modules/align_trim_1.nf
process ALIGN_TRIM_1 {
	tag { sampleId }
  publishDir "${params.out_dir}/medaka", mode: 'copy'

  input:
    tuple val(sampleId), path(input_bam), val(item), val(scheme), val(version), path(bed), path(align_trim)

  output:
    path "${params.run_name}_${sampleId}.alignreport_1.txt",														emit: align_report
    tuple val(sampleId), path("${params.run_name}_${sampleId}.trimmed.rg.sorted.bam"),	emit: trimmed_bam
    path "${params.run_name}_${sampleId}.trimmed.rg.sorted.bam.bai",										emit: trimmed_bai

  maxForks params.queueSize

  script:
  """
  set -euo pipefail
  python "${align_trim}" --normalise ${params.medaka_normalise} "${bed}" --start \
    --report "${params.run_name}_${sampleId}.alignreport_1.txt" \
    < "${input_bam}" \
    2> "${params.run_name}_${sampleId}.alignreport.err" \
  | samtools sort -T "${sampleId}" -o "${params.run_name}_${sampleId}.trimmed.rg.sorted.bam"
  samtools index "${params.run_name}_${sampleId}.trimmed.rg.sorted.bam"
  """
}


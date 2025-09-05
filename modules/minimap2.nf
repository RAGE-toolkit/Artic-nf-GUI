// modules/minimap2.nf
process MINIMAP2 {

  tag { sampleId }

  publishDir "${params.out_dir}/medaka", mode: 'copy', overwrite: true

  input: 
		tuple val(sampleId), path(reads), val(item), val(scheme), val(version), path(reference)

  output:
    tuple val(sampleId), path("${params.run_name}_${sampleId}.sorted.bam"),     emit: sorted_bam
    path "${params.run_name}_${sampleId}.sorted.bam.bai", emit: sorted_bai

  script:
    """
    set -euo pipefail

    minimap2 -a -x map-ont -t ${params.threads} \
			"${reference}"	\
			${reads.join(" ")}	\
    | samtools view -bS -F 4 - \
    | samtools sort -o "${params.run_name}_${sampleId}.sorted.bam"

    samtools index "${params.run_name}_${sampleId}.sorted.bam"
    """
}

//"${reference}" "${reads}"

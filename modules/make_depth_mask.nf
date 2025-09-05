// modules/make_depth_mask.nf
process MAKE_DEPTH_MASK {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_vcf), path(input_bam), path(reference), path(make_depth_mask_script), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.coverage_mask.txt"), emit: coverage_mask

	// errorStrategy 'terminate'   // switch to 'ignore' only if you really want to continue on failure

	script:
		if( !make_depth_mask_script.exists() )
			exit 1, "MAKE_DEPTH_MASK: Script not found: ${make_depth_mask_script}"
		if( !reference.exists() )
			exit 1, "MAKE_DEPTH_MASK: Reference not found: ${reference}"
		if( !input_bam.exists() )
			exit 1, "MAKE_DEPTH_MASK: Input BAM not found: ${input_bam}"
		if( !input_vcf.exists() )
			log.warn "MAKE_DEPTH_MASK: PASS VCF not found (continuing without it): ${input_vcf}"

	"""
	set -euo pipefail
	[ -s "${input_bam}.bai" ] || samtools index "${input_bam}"

	python "${make_depth_mask_script}" \
		--depth ${params.mask_depth} \
		--store-rg-depths \
		"${reference}" \
		"${input_bam}" \
		"${params.run_name}_${sampleId}.coverage_mask.txt"
	"""
}


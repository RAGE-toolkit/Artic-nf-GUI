// modules/mask.nf
process MASK {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(coverage_mask), path(fail_vcf), path(reference), path(mask_script), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.preconsensus.fasta"), emit: preconsensus

	// errorStrategy 'terminate'   // switch to 'ignore' if you prefer to continue on failure

	script:

		if( !mask_script.exists() )      exit 1, "MASK: Script not found: ${mask_py}"
		if( !reference.exists() )        exit 1, "MASK: Reference not found: ${reference}"
		if( !coverage_mask.exists() )    exit 1, "MASK: Coverage mask not found: ${coverage_mask}"
		if( !fail_vcf.exists() )         log.warn "MASK: FAIL VCF not found; proceeding without it: ${fail_vcf}"

		"""
		set -euo pipefail

		python "${mask_script}" \
			"${reference}" \
			"${coverage_mask}" \
			"${fail_vcf}" \
			"${params.run_name}_${sampleId}.preconsensus.fasta"
		"""
}


// modules/concat_for_muscle.nf
process CONCAT_FOR_MUSCLE {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(consensus_fa), path(reference), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${sampleId}.muscle.in.fasta"), emit: muscle_fa

	// errorStrategy 'terminate'   // change to 'ignore' if you want to continue on failure

	script:
		if( !consensus_fa.exists() )
			exit 1, "CONCAT_FOR_MUSCLE: Consensus FASTA not found: ${consensus_fa}"
		if( !reference.exists() )
			exit 1, "CONCAT_FOR_MUSCLE: Reference FASTA not found: ${reference}"

		"""
		set -euo pipefail
		cat "${consensus_fa}" "${reference}" > "${sampleId}.muscle.in.fasta"
    """
}


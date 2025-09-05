// modules/fasta_header.nf
process FASTA_HEADER {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(consensus_fa), path(fasta_header_script), val(item), val(scheme), val(version)

  output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.consensus.fasta"), emit: fasta

	// errorStrategy 'terminate'   // change to 'ignore' if you prefer to continue on failure

	script:
		//def header_script = file("${projectDir}/scripts/fasta_header.py")

		if( !consensus_fa.exists() )
			exit 1, "FASTA_HEADER: Input consensus FASTA not found: ${consensus_fa}"
		if( !fasta_header_script.exists() )
			exit 1, "FASTA_HEADER: Script not found: ${fasta_header_script}"

		"""
		set -euo pipefail

		#cp "${consensus_fa}" "${sampleId}.consensus.fasta"

		python "${fasta_header_script}" "${params.run_name}_${sampleId}.consensus.fasta" "${sampleId}"
    """
}


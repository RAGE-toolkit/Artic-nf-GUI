// modules/concat.nf
process CONCAT {

	tag { "concat_genome" }
	
	publishDir "${params.out_dir}/concatenate", mode: 'copy'

	input:
		path fasta_files
		path concat_script

	output:
		path "concat_genome.fasta", emit: genome_fa

	script:
		if( fasta_files == null || fasta_files.size() == 0 )
			exit 1, "CONCAT: No input FASTA files provided."

		"""
			set -euo pipefail
			python ${concat_script} -i ${fasta_files} -o concat_genome.fasta
    """
}
// 21       #cat ${fasta_files.join(' ')} > concat_genome.fasta


// modules/medaka_SNP-2.nf
process MEDAKA_SNP_2 {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_hdf), path(reference), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.2.vcf"), emit: vcf

	// errorStrategy 'terminate'   // prefer to fail loudly; switch to 'ignore' only if you must

	script:
		def reference = file("${params.primer_schema}/${scheme}/${version}/${scheme}.reference.fasta")

		if( !input_hdf.exists() )
			exit 1, "MEDAKA_SNP_2: Input HDF not found: ${input_hdf}"
		if( !reference.exists() )
			exit 1, "MEDAKA_SNP_2: Reference not found: ${reference}"

		"""
		set -euo pipefail

		medaka snp "${reference}" "${input_hdf}" "${params.run_name}_${sampleId}.2.vcf"
		"""
}


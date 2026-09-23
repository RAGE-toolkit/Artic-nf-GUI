// modules/medaka_SNP-2.nf
process MEDAKA_SNP_2 {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_hdf), path(reference), path(input_bam), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.2.vcf"), emit: vcf

	// errorStrategy 'terminate'   // prefer to fail loudly; switch to 'ignore' only if you must

	script:
		def reference = file("${params.primer_schema}/${scheme}/${version}/${scheme}.reference.fasta")

		if( !input_hdf.exists() )
			exit 1, "MEDAKA_SNP_2: Input HDF not found: ${input_hdf}"
		if( !reference.exists() )
			exit 1, "MEDAKA_SNP_2: Reference not found: ${reference}"
		if( !input_bam.exists() )
			exit 1, "MEDAKA_SNP_2: Input BAM not found: ${input_bam}"

		"""
		set -euo pipefail
		[ -s "${input_bam}.bai" ] || samtools index "${input_bam}"

		medaka snp "${reference}" "${input_hdf}" "${params.run_name}_${sampleId}.2.vcf"

		medaka tools annotate --pad 25 --RG 2 \
			"${params.run_name}_${sampleId}.2.vcf" \
			"${reference}" \
			"${input_bam}" \
			"${params.run_name}_${sampleId}.2.annotate.tmp.vcf"

		mv "${params.run_name}_${sampleId}.2.annotate.tmp.vcf" "${params.run_name}_${sampleId}.2.vcf"
		"""
}


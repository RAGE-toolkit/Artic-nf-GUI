// modules/longshot.nf
//longshot.nf

process LONGSHOT {

	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(input_vcf), path(input_bam), path(reference), val(item), val(scheme), val(version)

	output:
		tuple val(sampleId), path ("${params.run_name}_${sampleId}.longshot.merged.vcf"), emit: vcf
		//path "${sampleId}.potential.vcf.gz", emit: potential_vcf

	script:
		if( !input_vcf.exists() )  exit 1, "LONGSHOT: Input VCF not found: ${input_vcf}"
		if( !input_bam.exists() )  exit 1, "LONGSHOT: Input BAM not found: ${input_bam}"
		if( !reference.exists() )  exit 1, "LONGSHOT: Reference not found: ${reference}"

	"""
	set -euo pipefail

	if [ ! -f "${reference}.fai" ]; then
		echo "Indexing reference with samtools faidx"
		samtools faidx "${reference}"
	fi

	if [ ! -f "${input_bam}.bai" ] && [ ! -f "${input_bam}.bai" ]; then
		echo "Indexing BAM"
		samtools index "${input_bam}"
	fi

	longshot -P 0 -F -A --no_haps \
		--bam "${input_bam}" \
		--ref "${reference}" \
		--out "${params.run_name}_${sampleId}.longshot.merged.vcf"	\
	"""
}

//--potential_variants "${sampleId}.potential.vcf.gz"


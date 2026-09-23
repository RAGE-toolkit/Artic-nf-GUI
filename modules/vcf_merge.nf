// modules/vcf_merge.nf

process VCF_MERGE {
	tag { sampleId }

	publishDir "${params.out_dir}/medaka", mode: 'copy'

	input:
		tuple val(sampleId), path(vcf2), path(vcf1), path(vcf_unmatched), path(bed), path(vcf_merge_script), val(item), val(scheme), val(version)
	output:
		tuple val(sampleId), path("${params.run_name}_${sampleId}.merged.vcf.gz"),     emit: merged_vcf
		tuple val(sampleId), path ("${params.run_name}_${sampleId}.merged.vcf.gz.tbi"), emit: merged_tbi
		path "${params.run_name}_${sampleId}.primersitereport.txt", emit: primer_report

	script:
	"""
	set -euo pipefail

	# Sanity checks
	[ -s "${vcf1}" ] || { echo "VCF1 missing: ${vcf1}" >&2; exit 1; }
	[ -s "${vcf2}" ] || { echo "VCF2 missing: ${vcf2}" >&2; exit 1; }
	[ -s "${vcf_unmatched}" ] || { echo "Unmatched VCF missing: ${vcf_unmatched}" >&2; exit 1; }
	[ -s "${bed}"  ] || { echo "BED missing: ${bed}"  >&2; exit 1; }

	# vcf_merge.py needs primalbedtools, which only lives in the clair3env
	# conda env (see Dockerfile) - the system "python" does not have it.
	/opt/miniforge/envs/clair3env/bin/python3 "${vcf_merge_script}" \
		"${params.run_name}_${sampleId}" \
		"${bed}" \
		"2:${vcf2}" \
		"1:${vcf1}" \
		"unmatched:${vcf_unmatched}" \
		2> "${params.run_name}_${sampleId}.primersitereport.txt"

	bgzip -f "${params.run_name}_${sampleId}.merged.vcf"
	tabix -f -p vcf "${params.run_name}_${sampleId}.merged.vcf.gz"
	"""
}

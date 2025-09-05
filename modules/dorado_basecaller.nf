//dorado_basecaller.nf

process DORADO_BASECALLER {

	label "dorado_basecaller"

	publishDir "${params.out_dir}/dorado_basecaller", mode: 'copy', overwrite: true

	input:
		path basecaller_dir
		path model_dir
		path rawfile_dir

	output:
		path "calls.fastq", emit: fastq

	script:
	"""
	set -euo pipefail

	#command -v ${basecaller_dir} >&2 || { echo "ERROR: dorado not found: ${basecaller_dir}" >&2; exit 127; }
	#test -d "${model_dir}" || { echo "ERROR: model directory not found: ${model_dir}" >&2; exit 2; }

	"${basecaller_dir}" basecaller -r \
		"${model_dir}/${params.basecaller_config}" \
		-x "${params.run_mode}" \
		--emit-fastq \
		"${rawfile_dir}" \
		> "calls.fastq"
	"""
}


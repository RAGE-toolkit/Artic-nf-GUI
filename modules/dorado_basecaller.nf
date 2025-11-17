//dorado_basecaller.nf
process DORADO_BASECALLER {

        label "dorado_basecaller"

        publishDir "${params.out_dir}/dorado_basecaller", mode: 'copy', overwrite: true

        input:
                tuple path(rawfile_dir), val(basecaller_dir), path(basecaller_model_dir), val(basecaller_model), val(gpu)

        output:
                path "calls.fastq", emit: fastq

        script:
        """
        ${basecaller_dir} basecaller -r \
                "${basecaller_model_dir}/${basecaller_model}" \
                -x "${gpu}" \
                --emit-fastq \
                "${rawfile_dir}" \
                > "calls.fastq"
        """
}



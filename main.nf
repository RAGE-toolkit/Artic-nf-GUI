// main.nf
nextflow.enable.dsl=2

// -----------------------------------------------------------------------------
// Project-relative includes (robust under EPI2ME instances)
// -----------------------------------------------------------------------------
def MODULES = "${projectDir}/modules"

include { GUPPY_BASECALLER }   from "${MODULES}/guppy_basecaller.nf"
include { GUPPY_BARCODER }     from "${MODULES}/guppy_barcode.nf"
include { GUPPY_PLEX }         from "${MODULES}/guppy_plex.nf"

include { DORADO_BASECALLER }  from "${MODULES}/dorado_basecaller.nf"
include { DORADO_BARCODER }    from "${MODULES}/dorado_barcoder.nf"

include { PLEX_FQ_FILES }      from "${MODULES}/plex_fq_files.nf"
include { PLEX_DIRS }          from "${MODULES}/plex_dirs.nf"

include { MINIMAP2 }           from "${MODULES}/minimap2.nf"
include { ALIGN_TRIM }         from "${MODULES}/align_trim.nf"
include { SPLIT_UNMATCHED }    from "${MODULES}/split_unmatched.nf"
include { CLAIR3 }             from "${MODULES}/clair3-1.nf"
include { CLAIR3_2 }           from "${MODULES}/clair3-2.nf"
include { CLAIR3_UNMATCHED }   from "${MODULES}/clair3_unmatched.nf"
include { VCF_MERGE }          from "${MODULES}/vcf_merge.nf"
include { VCF_FILTER }         from "${MODULES}/vcf_filter.nf"
include { COMPRESS_AND_INDEX_VCF } from "${MODULES}/compress_and_index_vcf.nf"
include { MAKE_DEPTH_MASK }    from "${MODULES}/make_depth_mask.nf"
include { MASK }               from "${MODULES}/mask.nf"
include { BCFTOOLS_NORM }      from "${MODULES}/bcftools_norm.nf"
include { BCFTOOLS_CONSENSUS } from "${MODULES}/bcftools_consensus.nf"
include { FASTA_HEADER }       from "${MODULES}/fasta_header.nf"
include { CONCAT }             from "${MODULES}/concat.nf"
include { MAFFT }              from "${MODULES}/mafft.nf"
include { SUMMARY_STATS }      from "${MODULES}/summary_stats.nf"
include { REPORT }             from "${MODULES}/report.nf"

// -----------------------------------------------------------------------------
// Param defaults (CLI/GUI can override). No filesystem ops at compile time.
// -----------------------------------------------------------------------------
params.out_dir       = params.output_dir    ?: "${projectDir}/results"
params.fastq_dir     = params.fastq_dir     ?: 'raw_files/fastq'
params.rawfile_dir   = params.rawfile_dir   ?: 'raw_files'
params.rawfile_type  = params.rawfile_type  ?: 'fastq'              // 'fastq' | 'fast5_pod5'
params.basecaller    = params.basecaller    ?: 'Dorado'             // 'Dorado' | 'Guppy'
params.fq_extension  = params.fq_extension  ?: '.fastq'
params.threads       = (params.threads ?: 5) as int

// Common model/config params should already be set in nextflow.config (as you have)

// -----------------------------------------------------------------------------
// Sample sheet: require GUI/CLI to provide --sample_sheet (CSV with header)
// Expect columns: sampleId, barcode, schema, version
// -----------------------------------------------------------------------------
if( !params.sample_sheet ) {
  exit 1, "No sample sheet provided. Set --sample_sheet via EPI2ME GUI or CLI."
}
def sampleSheet = file(params.sample_sheet)
if( !sampleSheet.exists() ) {
  exit 1, "Sample sheet not found: ${sampleSheet}"
}

// Build CSV channel (light validation)
Channel
  .fromPath(sampleSheet, checkIfExists: true)
  .splitCsv(header: true, sep: ',')
  .map { row ->
      def sid     = (row.sampleId ?: '').toString().trim()
      def barcode = (row.barcode  ?: '').toString().trim()
      def scheme  = (row.schema   ?: '').toString().trim()
      def version = (row.version  ?: '').toString().trim()
      if( !sid || !barcode || !scheme || !version )
          throw new IllegalArgumentException("Sample sheet row missing fields: ${row}")
      tuple(sid, barcode, scheme, version)
  }
  .set { fq_channel }

// -----------------------------------------------------------------------------
// Logs (no mkdirs/prefetching here; processes handle their own outputs)
// -----------------------------------------------------------------------------
log.info ""
log.info "=== RAGE-toolkit/Artic-nf ==="
log.info "Output dir       : ${params.out_dir}"
log.info "Sample sheet     : ${sampleSheet}"
log.info "Run name         : ${params.run_name}"
log.info "Rawfile type     : ${params.rawfile_type}"
log.info "Basecaller       : ${params.basecaller}"
log.info "FASTQ dir        : ${params.fastq_dir}"
log.info "Rawfile dir      : ${params.rawfile_dir}"
log.info "Threads          : ${params.threads}"
log.info "QueueSize        : ${params.queueSize}"
log.info "Mask size        : ${params.mask_depth}"
log.info "Normalise        : ${params.normalise}"
log.info "Clair3 model     : ${params.model_path}"
log.info "Sequence length  : ${params.seq_len}"

if( params.rawfile_type == 'fast5_pod5' ) {
    log.info "---- Basecaller details ----"
    log.info "Basecaller config     : ${params.basecaller_config}"
    log.info "GPU mode              : ${params.gpu_mode}"
    if( params.basecaller_dir ) {
        log.info "Basecaller dir    : ${params.basecaller_dir}"
    }
    if( params.model_dir ) {
        log.info "Model dir         : ${params.model_dir}"
    }
}

log.info "=============================="

def dir_plex_script          = file("${projectDir}/scripts/plex.py", checkIfExists: true)
def vcf_merge_script	     = Channel.fromPath("${projectDir}/scripts/vcf_merge.py")
def vcf_filter_script        = Channel.fromPath("${projectDir}/scripts/vcf_filter.py")
def make_depth_mask_script   = Channel.fromPath("${projectDir}/scripts/make_depth_mask.py")
def mask_script              = Channel.fromPath("${projectDir}/scripts/mask.py")
def fasta_header_script      = Channel.fromPath("${projectDir}/scripts/fasta_header.py")
def summary_stats_script     = Channel.fromPath("${projectDir}/scripts/summary_stats.py")
def report_script            = Channel.fromPath("${projectDir}/scripts/report.py")
def concat_script            = Channel.fromPath("${projectDir}/scripts/consensus_combiner.py")


def medaka_dir 				= Channel.fromPath("${params.out_dir}/medaka")
def summary_stats_dir 		= Channel.fromPath("${params.out_dir}/summary_stats")
def rawfile_dir 			= file("${params.rawfile_dir}")
def basecaller_dir 			= Channel.fromPath("${params.basecaller_dir}")
def basecaller_model 		= Channel.fromPath("${params.model_dir}")

//================output directory========================

def ref_ch = fq_channel.map { sid, item, scheme, version ->
	def ref = file("${params.primer_schema}/${scheme}/${version}/${scheme}.reference.fasta")
	if( !ref.exists() ) throw new IllegalArgumentException("Missing reference for ${scheme}/${version} -> ${ref}")
	tuple(sid, ref)
	}

def bed_ch = fq_channel.map { sid, item, scheme, version ->
	def bed = file("${params.primer_schema}/${scheme}/${version}/${scheme}.scheme.bed")
	if( !bed.exists() ) throw new IllegalArgumentException("Missing bed file for ${scheme}/${version} -> ${bed}")
	tuple(sid, bed)
	}

plex_dirs_channel = fq_channel.map { sample_id, item, scheme, version ->
	tuple(rawfile_dir, dir_plex_script, sample_id, item, scheme, version)
	}


vcf_filter_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/vcf_filter.py", checkIfExists: true)
	tuple(sid, script)
	}

make_depth_mask_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/make_depth_mask.py", checkIfExists: true)
	tuple(sid, script)
	}

mask_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/mask.py", checkIfExists: true)
	tuple(sid, script)
	}

fasta_header_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/fasta_header.py", checkIfExists: true)
	tuple(sid, script)
	}

// -----------------------------------------------------------------------------
//  - MINIMAP2 and downstream take the prepared inputs plus fq_channel.
// -----------------------------------------------------------------------------
workflow {

  if( params.rawfile_type == 'fastq' ) {

		//running PLEX_DIRS here
		PLEX_DIRS(plex_dirs_channel)

		reads_ch = PLEX_DIRS.out.reads
    .map { file ->
        def base     = file.baseName
        def sampleId = base.replaceFirst(/_barcode.*/, "")
        tuple(sampleId, file)
    }
    .groupTuple()

		meta_ch = fq_channel.map { sid, item, scheme, version ->
    tuple(sid, item, scheme, version)
		}

		// MINIMAP channel here
		minimap_channel = reads_ch
    	.join(meta_ch)
    	.join(ref_ch)

		//**********running MINIMAP here
		MINIMAP2(minimap_channel)

		//ALIGN_TRIM channel
		align_trim_channel = MINIMAP2.out.sorted_bam
				.join(meta_ch)   // → [sid, bam, sid, item, scheme, version]
				.join(bed_ch)    // → [sid, bam, sid, item, scheme, version, sid, bed]
	}
  else {
    if( params.basecaller == 'Dorado' ) {
      //DORADO_BASECALLER(basecaller_dir, basecaller_model, rawfile_dir)
      //DORADO_BARCODER(fastq_file: DORADO_BASECALLER.out)
      //PLEX_FQ_FILES(DORADO_BARCODER.out, fq_channel)
      //MINIMAP2(input_dir: PLEX_FQ_FILES.out.collect(), fq_channel)
    }
    else {
      // Guppy path
      //GUPPY_BASECALLER(fast5_or_pod5_dir: params.rawfile_dir)
      //GUPPY_BARCODER(fastq_file: GUPPY_BASECALLER.out)
      //PLEX_DIRS(input_dir: GUPPY_BARCODER.out, fq_channel)
      //MINIMAP2(input_dir: PLEX_DIRS.out.collect(), fq_channel)
    }
  }

	//**********running ALIGN_TRIM
	ALIGN_TRIM(align_trim_channel)

	//SPLIT_UNMATCHED channel
	split_unmatched_channel = ALIGN_TRIM.out.primertrimmed_bam
		.join(meta_ch)

	//**********running SPLIT_UNMATCHED
	SPLIT_UNMATCHED(split_unmatched_channel)

	//CLAIR3 (pool 1) channel
	clair3_1_channel = ALIGN_TRIM.out.pool1_bam
		.join(ref_ch)
		.join(meta_ch)

	//**********running CLAIR3 (pool 1)
	CLAIR3(clair3_1_channel)

	//CLAIR3_2 (pool 2) channel
	clair3_2_channel = ALIGN_TRIM.out.pool2_bam
		.join(ref_ch)
		.join(meta_ch)

	//**********running CLAIR3_2 (pool 2)
	CLAIR3_2(clair3_2_channel)

	//CLAIR3_UNMATCHED channel
	clair3_unmatched_channel = SPLIT_UNMATCHED.out.unmatched_bam
		.join(ref_ch)
		.join(meta_ch)

	//**********running CLAIR3_UNMATCHED
	CLAIR3_UNMATCHED(clair3_unmatched_channel)

	//channel the vcf_merge script
	vcf_merge_scr = fq_channel.map { sid, item, scheme, version ->
	def script = file("${projectDir}/scripts/vcf_merge.py", checkIfExists: true)
		tuple(sid, script)
	}

	//VCF_MERGE channel
	vcf_merge_channel = CLAIR3_2.out.vcf
		.join(CLAIR3.out.vcf)
		.join(CLAIR3_UNMATCHED.out.vcf)
		.join(bed_ch)
		.join(vcf_merge_scr)
		.join(fq_channel)

	//**********running VCF_MERGE
	VCF_MERGE(vcf_merge_channel)

	//VCF_FILTER channel
	vcf_filter_channel = VCF_MERGE.out.merged_vcf
		.join(vcf_filter_scr)
		.join(fq_channel)

	//**********running VCF_FILTER
	VCF_FILTER(vcf_filter_channel)

	//COMPRESS_AND_INDEX_VCF channel
	compress_and_index_vcf_channel = VCF_FILTER.out.pass_vcf
		.join(fq_channel)

	//**********running COMPRESS_AND_INDEX_VCF
	COMPRESS_AND_INDEX_VCF(compress_and_index_vcf_channel)

	//MAKE_DEPTH_MASK channel
	make_depth_mask_channel = VCF_FILTER.out.pass_vcf
		.join(ALIGN_TRIM.out.primertrimmed_bam)
		.join(ref_ch)
		.join(make_depth_mask_scr)
		.join(fq_channel)

	//**********running MAKE_DEPTH_MASK
	MAKE_DEPTH_MASK(make_depth_mask_channel)

	//MASK channel
	mask_channel = MAKE_DEPTH_MASK.out.coverage_mask
		.join(VCF_FILTER.out.fail_vcf)
		.join(ref_ch)
		.join(mask_scr)
		.join(fq_channel)

	//**********running MASK
	MASK(mask_channel)

	//BCFTOOLS_NORM channel
	bcftools_norm_channel = MASK.out.preconsensus
		.join(COMPRESS_AND_INDEX_VCF.out.vcf_gz)
		.join(fq_channel)

	//**********running BCFTOOLS_NORM
	BCFTOOLS_NORM(bcftools_norm_channel)

	//BCFTOOLS_CONSENSUS channel
	bcftools_consensus_channel = MASK.out.preconsensus
		.join(BCFTOOLS_NORM.out.normalised_vcf)
		.join(BCFTOOLS_NORM.out.normalised_tbi)
		.join(MAKE_DEPTH_MASK.out.coverage_mask)
		.join(fq_channel)

	//**********running BCFTOOLS_CONSENSUS
	BCFTOOLS_CONSENSUS(bcftools_consensus_channel)


	//FASTA_HEADER channel
	fasta_header_channel = BCFTOOLS_CONSENSUS.out.consensus_fa
		.join(fasta_header_scr)
		.join(fq_channel)

	//**********running FASTA_HEADER
	FASTA_HEADER(fasta_header_channel)


	//CONCAT channel
	concat_channel = FASTA_HEADER.out.fasta
		.map { sid, fasta -> fasta }
		.collect()

	CONCAT(concat_channel, concat_script)

	//MAFFT channel
	mafft_channel = CONCAT.out.genome_fa

	//**********running MAFFT
	MAFFT(mafft_channel)

	//**********running SUMMARY_STATS
	SUMMARY_STATS(MAFFT.out.mafft_fa, medaka_dir, summary_stats_script)

	//**********running REPORT
	REPORT(SUMMARY_STATS.out.summary, medaka_dir, summary_stats_dir, report_script)
}

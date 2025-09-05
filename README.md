# 🧬 Artic-nf-GUI (RAGE-toolkit)

Artic-nf-GUI is a graphical user interface (GUI) for the [Artic-nf](https://github.com/RAGE-toolkit/Artic-nf) workflow, designed to integrate with Oxford Nanopore’s EPI2ME platform.

⚠️ Note: This project is currently under active development. At present, basecalling modules such as Dorado and Guppy are not supported within the GUI. Users should provide pre-basecalled FASTQ files (e.g., from a fastq_pass directory) as input to run the workflow.

The pipeline automatically selects the correct container image depending on your system architecture (`arm64/aarch64` or `x86_64`).  

---

## 📦 Requirements

- ⚠️ Note: This workflow has been developed and tested using **EPI2ME Desktop v5.2.5** on an **Apple M3 Pro system**. It is expected that other Apple Silicon (ARM64) platforms will also support successful execution.
A corresponding Docker image is available for Linux (Ubuntu), and while the workflow is likely to run on Linux systems, this has not yet been formally tested.

---

## 🚀 Quick start

### Download EPI2ME

Download the EPI2ME Desktop application here:
👉 https://epi2me.nanoporetech.com

### Import workflow in EPI2ME

Paste below Github link in the EPI2ME import
👉 https://github.com/RAGE-toolkit/Artic-nf-GUI

---

## Workflow parameters

Below are the main configurable parameters (default values shown):

| Parameter             | Default value                                 | Description                                              |
|-----------------------|-----------------------------------------------|----------------------------------------------------------|
| `--sample_sheet`      | `meta_data/meta_sheet.csv`                    | CSV file with sample metadata (**required**)             |
| `--output_dir`        | `results`                                     | Output directory                                         |
| `--run_name`          | `test_run`                                    | Run name (prefix for outputs)                            |
| `--rawfile_type`      | `fastq`                                       | Input type: `fastq` or `fast5_pod5`                      |
| `--rawfile_dir`       | `test_data/fastq_pass`                        | Input directory for raw files                            |
| `--fastq_dir`         | `raw_files/fastq`                             | Location for processed FASTQ files                       |
| `--fq_extension`      | `.fastq`                                      | FASTQ file extension                                     |
| `--primer_schema`     | `meta_data/primer-schemes`                    | Path to primer schemes                                   |
| `--kit_name`          | `EXP-NBD196`                                  | Nanopore kit name (e.g., `EXP-NBD104`, `EXP-NBD114`)     |
| `--threads`           | `5`                                           | CPUs for processes (minimap2, medaka, etc.)              |
| `--basecaller`        | `Dorado`                                      | Basecaller (`Dorado` or `Guppy`)                         |
| `--basecaller_dir`    | `null`                                        | Path to basecaller (only if not using container)         |
| `--model_dir`         | `null`                                        | Path to models (Dorado only)                             |
| `--basecaller_config` | `dna_r10.4.1_e8.2_400bps_fast@v4.2.0`         | Basecaller model config                                  |
| `--run_mode`          | `cuda:all`                                    | GPU run mode                                             |
| `--basecaller_threads`| `5`                                           | Guppy threads (ignored for Dorado)                       |
| `--medaka_model`      | `r941_min_fast_g303`                          | Medaka model                                             |
| `--medaka_normalise`  | `200`                                         | Medaka normalisation parameter                           |
| `--mask_depth`        | `20`                                          | Minimum depth for masking                                |
| `--seq_len`           | `350`                                         | Expected amplicon sequence length                        |


## Parameter description

The workflow is described by a [JSON schema](nextflow_schema.json), which defines all available parameters, their defaults, and validation rules.  
This schema is used to generate the **GUI configuration in EPI2ME Desktop** and ensures consistent parameter handling.

---

### 🔹 Input

- **meta_file** (`string`, *file-path*)  
  Path to the metadata file (CSV).  
  - Must be in CSV format  
  - See [📑 Sample sheet](meta_data/meta_sheet.csv) for details  

- **rawfile_dir** (`string`, *directory-path*)  
  Directory containing raw input files.  
  - Supports `FASTQ`, `FAST5`, or `POD5`  

- **rawfile_type** (`string`)  
  Type of raw input files.  
  - Options: `fastq`, `fast5_pod5`  
  - Default: `fastq`  

- **primer_schema** (`string`, *directory-path*)  
  Path to the primer scheme directory.  
  - Example: `/Documents/GitHub/Artic-nf/meta_data/primer_scheme`  

- **basecaller** (`string`)  
  Basecaller to run.  
  - Options: `Dorado`, `Guppy`  
  - Default: `Dorado`  

- **kit_name** (`string`)  
  Default: `EXP-NBD196`  
  - Defines barcode kit to be used  
  - Barcode lists can be obtained from Guppy/Dorado repositories  

- **run_name** (`string`)  
  Default: `test_run`  
  - Label for the sequencing run  
  - Used as a prefix for output files  

---

### 🔹 Output Options

- **output_dir** (`string`, *directory-path*)  
  Directory where all workflow results are stored.  
  - Default: `results`  

---

### 🔹 Basecalling Options

- **basecaller_dir** (`string`, *file-path*)  
  Path to the basecaller executable (optional).  

- **model_dir** (`string`, *directory-path*)  
  Path to basecalling models (required only for Dorado).  

- **basecaller_config** (`string`)  
  Default: `dna_r10.4.1_e8.2_400bps_fast@v4.2.0`  
  - Model config to use for basecalling.  

- **basecaller_threads** (`integer`)  
  Default: `5`  
  - Number of CPU threads for Guppy (ignored for Dorado).  

- **gpu_mode** (`string`)  
  Default: `cuda:all`  
  - GPU device specification for basecalling.  
  - Options: `cuda:0`, `cuda:all`, or `none`.  

---

### 🔹 Advanced Options

- **fastq_dir** (`string`)  
  Default: `raw_files/fastq`  
  - Directory for guppyplexed FASTQ files.  

- **fq_extension** (`string`)  
  Default: `.fastq`  
  - FASTQ file extension.  

- **seq_len** (`integer`)  
  Default: `350`  
  - Expected sequence length for plexing (primer dependent).  

- **medaka_model** (`string`)  
  Default: `r941_min_fast_g303`  
  - Medaka model for consensus polishing.  

- **medaka_normalise** (`integer`)  
  Default: `200`  
  - Target depth for Medaka normalization.  
  - Subsamples reads to reduce runtime and memory usage.  

- **mask_depth** (`integer`)  
  Default: `20`  
  - Minimum depth threshold for consensus masking.  

---

### 🔹 Miscellaneous Options

- **threads** (`integer`)  
  Default: `5`  
  - Number of CPU threads for multiprocess-enabled steps.  
  - ⚠️ Note: Minimap2 memory usage scales with thread count.  

- **queueSize** (`integer`)  
  Default: `5`  
  - Maximum number of parallel tasks allowed by the executor.  
  - Set to `1` for serial debugging.  

---

### 💻 Resources

- **Recommended**: 5 CPUs, 10 GB RAM  
- **Minimum**: 4 CPUs, 8 GB RAM  
- **Runtime**: ~1 minute per sample (depends on read count, reference length, and compute power)  
- **Architecture Support**: ✅ ARM64 (Apple Silicon, aarch64)  

---

📌 For the full schema, see [`nextflow_schema.json`](nextflow_schema.json).



# 🧬 Artic-nf (RAGE-toolkit)

Artic-nf-GUI is a graphical user interface (GUI) for the [Artic-nf]() workflow, designed to integrate with Oxford Nanopore’s EPI2ME platform.

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


### 📑 Sample sheet

The `sample_sheet` parameter should point to a CSV file containing the list of samples and their corresponding metadata such as **barcode**, **primer scheme**, and **version**.  

The file must include the following headers:


### Example

```csv
sampleId,barcode,schema,version
sampleA,barcode01,RABV,1
sampleB,barcode02,RABV,1

sampleId → Unique identifier for the sample
barcode → Barcode name used for demultiplexing (e.g., barcode01)
schema → Primer scheme (must match a scheme available in meta_data/primer-schemes/)
version → Version of the primer scheme
📌 For a complete example, see the provided file: meta_data/meta_sheet.csv

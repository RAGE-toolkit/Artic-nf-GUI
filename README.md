# 🧬 Artic-nf (RAGE-toolkit)

Artic-nf-GUI is a graphical user interface (GUI) for the Artic-nf workflow, designed to integrate with Oxford Nanopore’s EPI2ME platform.
⚠️ Note: This project is currently under active development. At present, basecalling modules such as Dorado and Guppy are not supported within the GUI. Users should provide pre-basecalled FASTQ files (e.g., from a fastq_pass directory) as input to run the workflow.

The pipeline automatically selects the correct container image depending on your system architecture (`arm64/aarch64` or `x86_64`).  

---

## 📦 Requirements

- [Nextflow ≥ 23.10.0](https://www.nextflow.io/)  
- Either:
  - **Docker** (recommended)  
  - or **Conda** (requires `mamba`/`conda`)  

---

## 🚀 Quick start

### Clone the repository

```bash
git clone https://github.com/sandeepkasaragod/test_gui.git
cd test_gui


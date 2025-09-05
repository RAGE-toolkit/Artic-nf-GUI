# 🧬 Artic-nf (RAGE-toolkit)

This pipeline performs **basecalling, demultiplexing, alignment, variant calling, consensus generation, and multiple sequence alignment** for viral genomes from **Nanopore sequencing data**.  

It supports both **Dorado** and **Guppy** basecallers, and can run either with **Docker** or **Conda** environments.  
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


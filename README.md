# 🧬 Artic-nf (RAGE-toolkit)

Artic-nf-GUI is a graphical user interface (GUI) for the Artic-nf workflow, designed to integrate with Oxford Nanopore’s EPI2ME platform.

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


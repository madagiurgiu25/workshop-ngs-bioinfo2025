# NGS Data Analysis Tutorial

In this series of tutorials, you will learn step-by-step how to analyze NGS (Next-Generation Sequencing) data, including QC, alignment, variant calling and variant annotation. 

---

## Table of Contents

### Tutorials
1. [Setup and installation](tutorials/tutorial_install.md)  
   - Download and organize raw FASTQ data
   - Download the human reference genome
   - Install necessary bioinformatics tools

2. [From FASTQ to BAM: Preprocessing NGS Data](tutorials/tutorial_from_fastq_to_bam.md)  
   Learn how to:
   - Perform initial quality control (FastQC & MultiQC)
   - Trim reads using Trim Galore!
   - Index the reference genome and prepare mapping
   - Align reads using BWA
   - Sort, index, and generate BAM files
   - Compute coverage tracks
   - Deduplicate reads and keep properly paired reads
   - Run final QC on cleaned BAM files

---

## 📂 Repository Structure

```text
workshop-ngs-bioinfo2025/
├── tutorials/
│   ├── tutorial_install.md
│   └── tutorial_from_fastq_to_bam.md
├── data/
│   ├── raw_fastq/
│   ├── trimmed_fastq/
│   ├── index/
│   ├── mapping/
│   ├── dedup/
│   └── ppaired/
└── README.md

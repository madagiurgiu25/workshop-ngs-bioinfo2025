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

3. [From BAM to VCF](tutorials/tutorial_from_bam_to_vcf.md)  
   Learn how to:
   - Perform variant calling
   - Perform variant filtering
   - Perform QC: variants stats and compute on-target coverage
---

## Repository Structure

```text
workshop-ngs-bioinfo2025/
├── tutorials/
│   ├── tutorial_install.md
│   ├── tutorial_from_fastq_to_bam.md
│   ├── tutorial_from_bam_to_vcf.md
│   └── tutorial_workflow_automation.md
├── data/
│   ├── DMD_exons_all.bed
│   ├── raw_fastq/     
│   ├── trimmed_fastq/
│   ├── index/
│   ├── mapping/
│   ├── dedup/
│   └── ppaired/
    ....
│   ├── Snakefile
│   ├── config.yaml
└── README.md
```


The content of the subfolders in `data` will be generated during the tutorial.

## License 

This project is licensed under the MIT License.
You are free to use, modify, and distribute this software, provided that the original license and copyright notice are included. See the full license text in the [LICENSE](LICENSE) file.


## Contributors

This repository was part of the **Genomic variants in Mendelian disease** organized by Prof. Horia Banciu at UBB, Cluj, as part of Bioinformatics Master curricula. Trainers: Dr. Anna Benet Pagès, Dr. Andreas Laner, Mădălina Giurgiu-Kraljič.

The code and text was ChatGPT assisted.



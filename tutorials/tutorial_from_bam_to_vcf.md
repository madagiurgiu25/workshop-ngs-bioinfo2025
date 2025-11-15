# GATK Germline Variant Calling & Filtering Analysis Pipeline

This document describes a complete workflow for variant calling, filtering, comparison, and coverage analysis using GATK, bcftools, and bedtools.

Here’s your pipeline with a **clean Markdown Table of Contents (TOC)** that works in GitHub, GitLab, Obsidian, VS Code, and Jupyter notebooks.
You can paste this directly into a `.md` file — all section links will work automatically.

---

## Table of Contents

* [1. Prepare GATK Metadata](#1-prepare-gatk-metadata)
* [2. Variant Calling with GATK HaplotypeCaller](#2-variant-calling-with-gatk-haplotypecaller)
* [3. Variant Filtering](#3-variant-filtering)
  * [3.1 Apply Quality Filters (VariantFiltration)](#31-apply-quality-filters-variantfiltration)
  * [3.2 Select Variants (SelectVariants)](#32-select-variants-selectvariants)
* [4. Variant Comparison and Statistics](#4-variant-comparison-and-statistics)
  * [4.1 Merge Variants](#41-merge-variants)
  * [4.2 Generate Variant Statistics](#42-generate-variant-statistics)
* [5. Compute On-Target Coverage](#5-compute-on-target-coverage)
  * [5.1 Extract Exons for a Gene](#51-extract-exons-for-a-gene)
  * [5.2 Compute Mean Exon Coverage](#52-compute-mean-exon-coverage)

---


## 1. Prepare Gatk metadata

```bash
ref=index/GRCh38.primary_assembly.genome.cut
gatk CreateSequenceDictionary -R $ref.fa -O $ref.dict
```


## 2. Variant calling using Gatk HaplotypeCaller

In this step we will call all possible variations relative to the reference genome.

```bash
mkdir -p variants
ref=index/GRCh38.primary_assembly.genome.cut.fa
sample=HG00171_female_ERR034564

gatk HaplotypeCaller \
  -R $ref  \
  -I ppaired/${sample}.proper.sorted.bam \
  -O variants/${sample}.vcf.gz \
  -A AlleleFraction
``` 

## 3. Variant filtering

### 3.1  Apply Quality Filters (VariantFiltration)

In this step we label variants which do not meet different quality criterias like:

- Insufficient coverage: depth of coverage `DP` < x (recommendation x=20)
- Strand bias, meaning the variant is observed only on forward or reverse strand; here two metrics are available, i.e. `FS > x` and `SOR > y` (recommendations are x=60, y=3.0)
- Tail distance bias, meaning the variant is observed only towards the read ends; here the metric is `ReadPosRankSum < x` (recommendation is x=-8.0)
- Low quality of the variant using `QUAL < x`
- Calls for with the AF is low and does not reflect a germline variant (`AF < 0.4 ` or `(AD[1]/(AD[0]+AD[1]) < 0.4)`)


```bash
ref=index/GRCh38.primary_assembly.genome.cut.fa
sample=HG00171_female_ERR034564
gatk VariantFiltration \
  -R $ref \
  -V variants/${sample}.vcf.gz \
  -O variants/${sample}_flagged.vcf.gz \
  --filter-name "LowDP" --filter-expression "DP < 20" \
  --filter-name "StrandBiasFS" --filter-expression "FS > 60.0" \
  --filter-name "StrandBiasSOR" --filter-expression "SOR > 3.0" \
  --filter-name "TailDistBias" --filter-expression "ReadPosRankSum < -8.0" \
  --filter-name "NonGermline" --filter-expression "(AD[1]/(AD[0]+AD[1]) < 0.4)" \
  --filter-name "LowQual" --filter-expression  "QUAL < 10"
```

### 3.2 Select Variants (SelectVariants)

This step will remove all variants having a tag other than `PASS`.

```bash
sample=HG00171_female_ERR034564

gatk SelectVariants -V variants/${sample}_flagged.vcf.gz --exclude-filtered -O variants/${sample}_filtered.vcf.gz
```

If the user would like to keep `PASS` and some filter `NonGermline`, one could run:

```bash
sample=HG00171_female_ERR034564

gatk SelectVariants \
-V variants/${sample}_flagged.vcf.gz 
--select-expression "vc.getFilter() == 'PASS' || vc.getFilter() == 'NonGermline'"\
-O variants/${sample}_filtered_keep_PASS_NonGermline.vcf.gz
```

Repeat Steps 2 and 3 for the `HG00152_male_SRR769545`.

## 4. Variant Comparison and Statistics

### 4.1 Merge variants

Using `bcftools merge` use the :

```bash
bcftools merge variants/HG00171_female_ERR034564_filtered.vcf.gz variants/HG00152_male_SRR769545_filtered.vcf.gz -o variants/merged_filtered.vcf.gz -O z
```

### 4.2 Generate Variant Statistics

Using `bcftools` create summary statistics of the SNVs and INDELs:

```bash
mkdir -p stats
sample=HG00171_female_ERR034564
bcftools stats variants/${sample}_filtered.vcf.gz > stats/${sample}.stats
plot-vcfstats stats/${sample}.stats -p stats_plots/
```


```bash
# compare two samples
bcftools isec -p isec_out variants/HG00171_female_ERR034564_filtered.vcf.gz variants/HG00152_male_SRR769545_filtered.vcf.gz

# create index
bgzip -c isec_out/0001.vcf > isec_out/0001.vcf.gz
bgzip -c isec_out/0002.vcf > isec_out/0002.vcf.gz
bgzip -c isec_out/0003.vcf > isec_out/0003.vcf.gz

tabix -p vcf isec_out/0001.vcf.gz
tabix -p vcf isec_out/0002.vcf.gz
tabix -p vcf isec_out/0003.vcf.gz
```

## 5. Compute on-target coverage

### 5.0 Required packages

Make sure you have installed:
- awk
- gawk
- bedtools

```
# install MacOS
brew install awk gawk bedtools

# install Linux
apt-get install awk gawk bedtools
```

### 5.1 Extract Exons for a Gene

To compute a BED file containing the union exons of the each target gene one can use the gene annotation from GENCODE. The code below assumes a basic familiarity with `awk` and `regex`.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.basic.annotation.gtf.gz
gunzip gencode.v49.primary_assembly.basic.annotation.gtf.gz
```

These are stored in [GTF](https://www.ensembl.org/info/website/upload/gff.html) format looking like:

```
chrX  HAVANA  gene  31097677  33339609  . - . gene_id "ENSG00000198947.18"; gene_type "protein_coding" gene_name "DMD";
chrX  HAVANA  exon  31177932  31178178  . - . gene_id "ENSG00000198947.18"; transcript_id "ENST00000679437.1"; gene_type "protein_coding"; gene_name "DMD";
chrX  HAVANA  exon  31173539  31173604  . - . gene_id "ENSG00000198947.18"; transcript_id "ENST00000679437.1"; gene_type "protein_coding"; gene_name "DMD";
...
```

Extract all the exons (including coding and non-coding exons) for DMD gene:

```bash

gene=DMD
gtf=gencode.v49.primary_assembly.basic.annotation.gtf

# Extract all lines in the file which match exactly "DMD" and are exons
# MacOS
grep -e "\""$gene"\"" $gtf | awk '$3=="exon"' > ${gene}_exons.gtf

# Linux
grep -P "\""$gene"\"" $gtf | awk '$3=="exon"' > ${gene}_exons.gtf
```

The  `${gene}_exons.gtf` contains the following information:

```bash
head -1 ${gene}_exons.gtf
chrX	HAVANA	exon	31177932	31178178	.	-	.	gene_id "ENSG00000198947.18"; transcript_id "ENST00000679437.1"; gene_type "protein_coding"; gene_name "DMD"; transcript_type "protein_coding"; transcript_name "DMD-231"; exon_number 1; exon_id "ENSE00003910833.1"; level 2; protein_id "ENSP00000506629.1"; hgnc_id "HGNC:2928"; tag "CAGE_supported_TSS"; tag "RNA_Seq_supported_only"; tag "basic"; tag "GENCODE_Primary"; havana_gene "OTTHUMG00000021336.7";
```

Next, extract per line the following information:

- Gene name using the regex `gene_name "([^"]+)";` - ATTENTION in this case `^` means negation, and `[^"]` means match everything which is not a quote `"`; `^` can also mean start of the line if it is not inside a `[ ]`
- Exon number using the regex `exon_number ([0-9]+);` - means match only numbers

```bash
# Convert to BED format
gawk -v gene="$gene" -F"\t" '
BEGIN{OFS="\t"}
{
    match($0,/gene_name "([^"]+)";/,g)
    match($0,/exon_number ([0-9]+);/,e)
    print $1,$4-1,$5,"exon_" e[1],".",$7,g[1]
}' "${gene}_exons.gtf" > "${gene}_exons.bed"
```

And output a BED format: `chr start end exon_number score strand gene`, which is the equivalent in your file of `print $1,$4-1,$5,"exon_"e[1],".",$7,g[1];`.

```bash
head ${gene}_exons.bed
# chrX	31177931	31178178	exon_1	.	-	DMD
# chrX	31173538	31173604	exon_2	.	-	DMD
# chrX	31172347	31172413	exon_3	.	-	DMD
```

Every gene has multiple isoforms. To compute the on-target coverage we are interested to compute the union of non-overlapping regions across all isoforms. For this we use `bedtools merge`. Upfront we need to sort all the coordinates. 

```bash
bedtools sort -i ${gene}_exons.bed | bedtools merge -i - -s -c 4,5,6,7 -o collapse,distinct,distinct,distinct > ${gene}_exons_all.bed
```

This will generate the following file:

```bash
head ${gene}_exons_all.bed
# chrX	31119218	31121930	exon_17,exon_8,exon_9,exon_51,exon_79,exon_25,exon_79,exon_34,exon_36,exon_32,exon_35,exon_31,exon_17,exon_14,exon_16,exon_15,exon_8,exon_33,exon_17,exon_18,exon_35,exon_16,exon_13	.	-	DMD
# chrX	31126641	31126673	exon_34,exon_16,exon_15,exon_32,exon_31,exon_78,exon_13,exon_16,exon_78,exon_50,exon_35,exon_24,exon_17	.	-	DMD
# chrX	31126885	31127186	exon_16	.	-	DMD
# chrX	31134101	31134194	exon_49,exon_14,exon_30,exon_7,exon_16,exon_12,exon_34,exon_16,exon_15,exon_77,exon_15,exon_33,exon_34,exon_77,exon_12,exon_15,exon_30,exon_15,exon_7,exon_23,exon_8,exon_33,exon_14	.	-	DMD
```

Finally we renumber column 4, i.e. `exon_1`, `exon_2`, ....

```bash
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"exon_"NR,$5,$6,$7}' ${gene}_exons_all.bed > tmp && mv tmp ${gene}_exons_all.bed

head ${gene}_exons_all.bed
# chrX	31119218	31121930	exon_1	.	-	DMD
# chrX	31126641	31126673	exon_2	.	-	DMD
# chrX	31126885	31127186	exon_3	.	-	DMD
# chrX	31134101	31134194	exon_4	.	-	DMD
```

Attention: Exons may be coding or non-coding, and when overlapping regions are merged, the exon numbering will no longer correspond to the exon numbers of individual isoforms.

### 5.2 Compute Mean Exon Coverage

We compute the on-target mean coverage of `DMD` gene for `HG00171_female_ERR034564` and `HG00152_male_SRR769545`. The mean coverage is compute by averaging the read-depth across all bases covering an exon.

```bash
# assumes you are in workshop-ngs-bioinfo2025/data
mkdir -p stats
gene=DMD

sample=HG00171_female_ERR034564
bedtools coverage -a index/${gene}_exons_all.bed -b ppaired/${sample}.proper.sorted.bam -mean > stats/${sample}_${gene}_exon_mean_coverage.bed

sample=HG00152_male_SRR769545
bedtools coverage -a index/${gene}_exons_all.bed -b ppaired/${sample}.proper.sorted.bam -mean > stats/${sample}_${gene}_exon_mean_coverage.bed
```

Output will contain an additional column wit the on-target mean coverage:

```bash
head stats/HG00171_female_ERR034564_DMD_exon_mean_coverage.bed 
# chrX	31097676	31098183	exon_1	.	-	DMD	0.9644970
# chrX	31119218	31121930	exon_2	.	-	DMD	17.6165199
# chrX	31126641	31126797	exon_3	.	-	DMD	181.5384674
# chrX	31126885	31127186	exon_4	.	-	DMD	9.0465117
# chrX	31129926	31134689	exon_5	.	-	DMD	28.1490650

head stats/HG00152_male_SRR769545_DMD_exon_mean_coverage.bed
# chrX	31097676	31098183	exon_1	.	-	DMD	0.3846154
# chrX	31119218	31121930	exon_2	.	-	DMD	1.1530236
# chrX	31126641	31126797	exon_3	.	-	DMD	32.5897446
# chrX	31126885	31127186	exon_4	.	-	DMD	0.0000000
# chrX	31129926	31134689	exon_5	.	-	DMD	1.3722444
```

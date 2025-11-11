

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

### 3.1 Set filters using Gatk VariantFiltration

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


### 3.2 Select variants using Gatk SelectVariants

This step will remove all variants having a tag other than `PASS`.

```bash
sample=HG00171_female_ERR034564

gatk SelectVariants -V variants/${sample}_flagged.vcf.gz --exclude-filtered -O variants/${sample}_filtered.vcf.gz
```

Repeat Steps 2 and 3 for the `HG00152_male_SRR769545`.

## 4. Variant comparison and stats

### 4.1 Merge variants

Using `bcftools merge` use the :

```bash
bcftools merge variants/HG00171_female_ERR034564_filtered.vcf.gz variants/HG00152_male_SRR769545_filtered.vcf.gz -o variants/merged_filtered.vcf.gz -O z
```

### 4.2 Stats

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

For a gene first extract all the exons:

```bash
# Extract all exons
gene=DMD
grep -w $gene anno.gtf | awk '$3=="exon"' > ${gene}_exons.gtf

# Convert to BED format
awk -v gene=$gene -F"\t" 'BEGIN{OFS="\t"} $3=="exon" && match($0,/gene_name "([^"]+)"/,g) && g[1]==gene {match($0,/gene_name "(\w+)";/,g); match($0,/exon_number ([0-9]+);/,e); print $1,$4-1,$5,"exon_"e[1],".",$7,g[1];}' ${gene}_exons.gtf > ${gene}_exons.bed

bedtools sort -i ${gene}_exons.bed | bedtools merge -i - -s -c 4,5,6,7 -o collapse,distinct,distinct,distinct > ${gene}_exons_all.bed

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"exon_"NR,$5,$6,$7}' ${gene}_exons_all.bed > tmp && mv tmp ${gene}_exons_all.bed
```

Compute mean coverage

```bash
gene=DMD

sample=HG00171_female_ERR034564
bedtools coverage -a index/${gene}_exons_all.bed -b ppaired/${sample}.proper.sorted.bam -mean > stats/${sample}_${gene}_exon_mean_coverage.bed

sample=HG00152_male_SRR769545
bedtools coverage -a index/${gene}_exons_all.bed -b ppaired/${sample}.proper.sorted.bam -mean > stats/${sample}_${gene}_exon_mean_coverage.bed
```


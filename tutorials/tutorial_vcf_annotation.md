# VCF Annotation

## 1. Installation

```bash
mamba create -n vep -c conda-forge -c bioconda ensembl-vep
```

## 2. Prepare Custome Genome

```
vep_install -a cf \
            -s GRCh38.primary_assembly.genome.cut \
            -y GRCh38.primary_assembly.genome.cut \
            -f index/GRCh38.primary_assembly.genome.cut.fa \
            --gtf index/anno.gtf \
            -c index/vep_custom_cache
```

```bash
docker run -t -i -v $PWD/variants:/mnt ensemblorg/ensembl-vep vep -i /mnt/HG00171_female_ERR034564_filtered.vcf.gz -o /mnt/HG00171_female_ERR034564_filtered.annotated.vcf.gz --vcf --everything --database --assembly GRCh38
```

```
# on cluster download human genome
/data/cephfs-1/work/groups/henssen/users/giurgium_c/miniforge3/envs/vep/share/ensembl-vep-115.2-1/vep_install  -a c -s homo_sapiens -y GRCh38 --AUTO c --CACHEDIR $PWD/.vep
```

```
./vep --af --af_gnomade --appris --biotype --buffer_size 500 --canonical --check_existing --check_frequency --custom file=https://ftp.ensembl.org[path_to]/gencode.v48.promoter_windows_sorted.gff3.gz,format=gff,short_name=GENCODE_promoter,type=overlap,gff_type=gencode_promoter --distance 5000 --freq_filter exclude --freq_freq 0.01 --freq_gt_lt gt --freq_pop 1kg_all --gencode_primary --hgvs --mane --numbers --pick --plugin SpliceAI,snv=[path_to]/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz,indel=[path_to]/spliceai_scores.masked.indel.hg38.vcf.gz,snv_ensembl=[path_to]/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz --plugin AlphaMissense,file=[path_to]/AlphaMissense_hg38.tsv.gz --plugin EVE,file=[path_to]/eve_merged.vcf.gz --plugin REVEL,file=[path_to]/new_tabbed_revel_grch38.tsv.gz --plugin GO,[path_to]/ --plugin NMD --plugin Enformer,file=[path_to]/enformer_grch38.vcf.gz --plugin CADD,snv=[path_to]/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz,indels=[path_to]/CADD_GRCh38_1.7_InDels.tsv.gz --polyphen b --pubmed --regulatory --show_ref_allele --sift b --species homo_sapiens --symbol --transcript_version --tsl --uploaded_allele --cache --input_file [input_data] --output_file [output_file]
```
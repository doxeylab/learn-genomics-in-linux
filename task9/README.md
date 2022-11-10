# Task9 - Analysis of human SNP variation from 1000 genomes data

In this lab, you will be analyzing available data from the 1000 genomes project - https://en.wikipedia.org/wiki/1000_Genomes_Project

You will 

### Requirements

#### Command-line tools
* Access to a linux-based OS running BASH
* Tabix
* vcftools


## Getting Started

* Login to your linux environment and create a new folder for your task7

```
mkdir task9  #creates folder
cd task9 #enters into folder
```

## Download the VCF file for your chromosome of interest

e.g., below we will download chromosome 12

```
ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
```

## Download the reference genome (optional, if needed)
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

gunzip human_g1k_v37.fasta.gz

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
```


## Use `tabix` to extract the region of interest from the chromosome 

e.g., suppose we are interested in the variants found across the 1000-bp region 49687909-49688909 of chromosome 12

```
tabix -fh ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 12:49687909-49687919 >region.vcf
```

## convert VCF file to tab-separated file
```
cat region.vcf | vcf-to-tab
```

* How many SNPs were detected?



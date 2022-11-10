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

## Retrieving the raw data

The `tabix` utility can be used to retrieve 1000 genomes data from the public FTP directory in the following way:

```
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 12:49687909-49689000
```


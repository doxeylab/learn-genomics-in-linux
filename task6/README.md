# Task6 - Resequencing: variant calling from NGS data

In this lab (based on this [tutorial](https://angus.readthedocs.io/en/2014/variant.html)) you will be exploring some resequencing data from Richard Lenski's famous E. coli long term evolution experiment.
More on Lenski's long term evolution experiment can be found in this [article](http://www.nature.com/nature/journal/v489/n7417/full/nature11514.html)

...and here: https://en.wikipedia.org/wiki/Richard_Lenski

You will be mapping the reads from a single population of E. coli at 38,000 generations that has evolved citrate-utilization capacity (Cit+). You will map the reads to the reference genome in order to identify SNPs that have occurred in this lineage (or its ancestors). You will focus on one SNP in particular that has created a "mutator" strain by disrupting the mutS DNA-repair gene.


### Requirements

* Access to a linux-based OS running BASH
* [bwa](http://bio-bwa.sourceforge.net/)
* [samtools](http://samtools.sourceforge.net/)
* [tablet](https://ics.hutton.ac.uk/tablet/download-tablet/)


## Getting Started

* Login to your linux environment and create a new folder for your task6

```
mkdir task6  #creates folder
cd task6 #enters into folder
```

## Retrieving the raw data

First, download the E. coli reference genome:

```
wget athyra.idyll.org/~t/REL606.fa.gz
gunzip REL606.fa.gz
```

Next, download the resequencing data. This is 229 MB so you may have to be patient.

```
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098038/SRR098038.fastq.gz
```


## Mapping your reads to the reference genome


Before we can map reads with `bwa` (note: `bowtie` is another option), we need to index the reference genome. This can be done with `bwa index`.

```
bwa index REL606.fa
```

Now, let's map the reads to the reference genome. This is also a fairly intensive step that may take a few minutes.

```
bwa aln REL606.fa SRR098038.fastq.gz > SRR098038.sai
```

Make a .SAM file which contains all information about where each read maps onto the reference genome

```
bwa samse REL606.fa SRR098038.sai SRR098038.fastq.gz > SRR098038.sam
```

Index the reference genome (again) so that `samtools` can work with it

```
samtools faidx REL606.fa
```

Convert .SAM file to .BAM file

```
samtools import REL606.fa.fai SRR098038.sam SRR098038.bam
```

Sort BAM file and index it

```
samtools sort SRR098038.bam SRR098038.sorted
samtools index SRR098038.sorted.bam
```

## Visualizing your mapped reads

Now, you can visualize your mapped reads and identify variants.

[<b>Download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) the following two files to your local machine. Tip: find the path to your file with `realpath yourFile.txt`

* SRR098038.sorted.bam
* REL606.fa

Then open them both in `tablet`.


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1 - Paste a screenshot of a portion of the mapped reads covering a few kilobases of DNA (pileup).

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2 - Use tablet to visually identify SNPs or other variants. Go to Adjust/Variants to turn up the contrast. Find a position that is likely a SNP and not a sequencing error. Paste a screenshot.


This lineage of E. coli has a mutation in the mutS gene. The mutation creates a premature stop codon. Can you find it!?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3 - Paste a screenshot of the region containing this mutation and what was the amino acid encoded by this codon before this mutation?

## Variant calling

Instead of identifying SNPs by eye, use `samtools` to call variants in an automated fashion.

```
samtools mpileup -uD -f REL606.fa SRR098038.sorted.bam | bcftools view -bvcg - > SRR098038.raw.bcf

#convert to vcf (human-readable variant call format). This file should contain all identified SNPs and other variants.
bcftools view -v -c -g SRR098038.raw.bcf > SRR098038.raw.vcf

```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q4 - How many total variants are present? Hint: use "grep -v"


---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) above.
Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.




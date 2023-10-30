# Task6 - Resequencing: variant calling from NGS data

In this lab (based on this [tutorial](https://angus.readthedocs.io/en/2014/variant.html)) you will be exploring some resequencing data from Richard Lenski's famous <i>E. coli</i> long term evolution experiment.
More on Lenski's long term evolution experiment can be found in this [article](http://www.nature.com/nature/journal/v489/n7417/full/nature11514.html), with a summary here: https://en.wikipedia.org/wiki/Richard_Lenski.

You will be mapping the reads from a single population of <i>E. coli</i> at 38,000 generations that has evolved citrate-utilization capacity (Cit+). You will map the reads to the reference genome in order to identify SNPs that have occurred in this lineage (or its ancestors). You will focus on one SNP in particular that has created a "mutator" strain by disrupting the <i>mutS</i> DNA-repair gene.


### Requirements

#### Command-line tools
* Access to a linux-based OS running BASH
* [bwa](http://bio-bwa.sourceforge.net/)
* [samtools](http://samtools.sourceforge.net/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)

#### Graphical tools

* You will also need to download and install either Tablet or IGV on your own machine.
* [tablet](https://ics.hutton.ac.uk/tablet/)
* [igv](http://software.broadinstitute.org/software/igv/)


## Getting Started

* Login to your linux environment and create a new folder for your task6.

```
mkdir task6  #creates folder
cd task6 #enters into folder
```

## Retrieving the raw data

First, download the <i>E. coli</i> reference genome:

```
wget https://raw.githubusercontent.com/doxeylab/learn-genomics-in-linux/master/task6/ecoli-rel606.fa.gz
gunzip ecoli-rel606.fa.gz
```

The resequencing data is located at `/data/SRR098038.fastq.gz`. This is a 229 MB file, so it was already downloaded for you using the following command:

```
#wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098038/SRR098038.fastq.gz
```


## Mapping your reads to the reference genome


Before we can map reads with `bwa` (note: `bowtie` is another popular option), we need to index the reference genome. This can be done with `bwa index`.

```
bwa index ecoli-rel606.fa
```

Now, let's map the reads to the reference genome. This is also a fairly intensive step that may take a few minutes.

```
bwa aln ecoli-rel606.fa /data/SRR098038.fastq.gz > SRR098038.sai
```

Make a .SAM file which contains all information about where each read maps onto the reference genome.

```
bwa samse ecoli-rel606.fa SRR098038.sai /data/SRR098038.fastq.gz > SRR098038.sam
```

Index the reference genome (again) so that `samtools` can work with it.

```
samtools faidx ecoli-rel606.fa
```

Convert the .SAM file to a .BAM file.

```
samtools view -S -b SRR098038.sam > SRR098038.bam
rm SRR098038.sam  # remove this large file
```

Sort the BAM file and index it.

```
samtools sort SRR098038.bam > SRR098038.sorted.bam
samtools index SRR098038.sorted.bam
```

## Viewing your BAM file

BAM files can be viewed with `igv` or with `tablet`. But let's take a quick look in the terminal.

```
samtools tview SRR098038.sorted.bam
```

Type `q` to exit when finished.

## Variant calling

Instead of identifying SNPs by eye, use `bcftools` to perform automated variant calling.

```
bcftools mpileup -f ecoli-rel606.fa SRR098038.sorted.bam | bcftools call -mv -Ob --ploidy 1 -o calls.bcf

#convert to vcf (human-readable variant call format). This file should contain all identified SNPs and other variants.
bcftools view calls.bcf > calls.vcf

```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1 - How many total variants are present? Hint: `grep` for a pattern found only in your variant lines.


## Locating a key SNP in Lenski's <i>E. coli</i> evolution experiment

This lineage of <i>E. coli</i> has a mutation in the <i>mutS</i> gene (protein sequence can be found [here](https://www.uniprot.org/uniprot/P23909.fasta)). This mutation creates a premature stop codon. Your task is to find this mutation within your sequencing data!

Find the region in the reference genome that encodes the <i>mutS</i> gene using `blast`. You may need to refer to earlier tasks to help you with this.

Now, extract the mapped area for this region from your .bam file. The command will be something like this:

```
samtools view SRR098038.sorted.bam "ecoli:START-END" > region.sam   # where START and END are position numbers
```


Download the following two files and open these files in `tablet` on your home machine.

- the region.sam file you created above
- the reference genome ([ecoli-rel606.fa](https://raw.githubusercontent.com/doxeylab/learn-genomics-in-linux/master/task6/ecoli-rel606.fa.gz))

Now, locate the region containing the <i>mutS</i> gene within `tablet`, and search for the premature stop codon variant.

Here is an [example read-pileup](https://raw.githubusercontent.com/doxeylab/learn-genomics-in-linux/master/task6/example-pileup.png) in tablet that highlights a variant position.


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2 - Paste a screenshot highlighting this mutation (you will need to zoom in) and show the amino acid translation<!--see [here](https://software.broadinstitute.org/software/igv/sequence_track_options)-->. 


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q3 - The premature stop codon mutation is from the codon “_ _ _ ” to “_ _ _ ”.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q4 - The amino acid encoded by this codon before this mutation is _____?


Once you are finished, please delete the files in your task 6 folder with the `rm` command.:

> :warning: **Caution:** Be careful as `rm` permanently deletes files!

```
cd task6
rm *
```



---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) above.
Submit your answers to questions 1-4 to the QUIZ on LEARN.




# Task3 - Genome Annotation

This task is a tutorial on genome annotation using `prokka` and other tools.

### Requirements

* Access to a linux-based OS running BASH
* Prokka

## Installation

If you do not already have access to a GUI running the graphical software listed above, please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

All software used are available for Mac/Windows/Linux.

---

## Getting Started

* Login to your linux environment as you did in task1.

* Create a new folder for your task3

```
mkdir task3  #creates folder
cd task3 #enters into folder
```

## Retrieving the raw data

* Copy the genome you assembled from task2

```
cp ../task2/abyss-assembly-contigs.fa . 
```

## Annotation of your genome from Task2

By marking the ORFs in your genome (given a min size threshold), you have essentially performed a simple gene finding algorithm. However, there are more advanced ways of gene-finding that take additional criteria into account.

A popular genome annotation tool for prokaryotic genomes is [`prokka`](https://github.com/tseemann/prokka).
`prokka` automates a series of genome annotation tools and is simple to run. It has been installed for you on the server.

* Run prokka using the following command

```
prokka abyss-assembly-contigs.fa
```

* Now, download the .gbk file that was produced and view it in Artemis

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1) Are the predicted gene locations consistent with your earlier ORF predictions from task 2? 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2) Why are there vertical black lines in the middle of predicted ORFs?

* Quit artemis and change your artemis 'Options' to better reflect the source of this genome. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3) Does this correct the issues above? Paste a screenshot.

## Annotation of an <i>E. coli</i> genome

Next, let's perform genome annotation on a larger scale.

* Download (or copy from your task1 folder) the E. coli K12 genome from task1 and run `prokka`

```
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/e-coli-k12-genome.fasta.gz
gunzip e-coli-k12-genome.fasta.gz
prokka /e-coli-k12-genome.fasta
```

Next, explore the files produced by `prokka`. Start with the .txt file.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q4) How many genes, rRNAs, tRNAs, and CRISPR loci were predicted? What is the size of the genome in Mb?

Prokka also annotates genes based on [COGs](https://www.ncbi.nlm.nih.gov/COG/) and also [E.C.](https://enzyme.expasy.org/) (enzyme commission) numbers. This information can be found in the .tbl file. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q5) How many genes were annotated with COGS? What proportion of total genes is this?

Column 6 of this .tsv file lists the COGs. To print out just column 6, do the following:

```
cut -f6 PROKKA_09182018.tsv
```

Next, let's sort this list using `sort` and pipe it to the `uniq` command.

```
cut -f6 PROKKA_09182018.tsv | sort | uniq
```

How many unique COGs were assigned? Just count the number of lines in the output

```
cut -f6 PROKKA_09182018.tsv | sort | uniq | wc -l
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q6) How many unique enzymatic activities (E.C. numbers) were assigned to the E. coli genome?

## Extracting specific regions of interest

Now, suppose you are interested in a non-coding (regulatory region) such as the promoter of the "trp operon". The trp operon promoter can be found directly upstream of the <b>trpE</b> gene.

Using tools such as `makeblastdb`, `blastdbcmd`, `grep`, `revseq` (which reverse complements a sequence if necessary), and the information in the .gff file, extract the 30-nucleotide long sequence upstream of the trpE gene.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q7) Paste this sequence into your assignment as well as the code you used to extract it.

### Advanced: Extracting the rRNAs predicted by barrnap

Here is a quick three-liner to extract the 16S rRNAs predicted by barrnap, and align them.

```
# get the regions of  from the .gff file
makeblastdb yourGenome.fna
cat *.gff | grep "barrnap" | awk '{ if ($7 == "-") {print $1" "$4"-"$5" minus"} else {print $1" "$4"-"$5" plus"} }' >rRNAs.txt
blastdbcmd -db yourGenome.fna -entry_batch list | muscle
```


## Recap: Annotating a novel genome of unknown source

And now for something a little more difficult.

Next, we will be giving you the raw reads for a sequencing project of an unknown organism (your only hint is that it is bacterial). Using any of the tools above: 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q7) Assemble the genome into a single contig, annotate the genes, and produce a circular plot of the genome using `dnaPlotter` in the `artemis` package. Include your source code.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q8) Based on 16S rRNA, what is the taxonomic source of this DNA? Include your source code.








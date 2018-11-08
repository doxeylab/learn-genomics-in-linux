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

Before we can map reads with `bwa` (note: `bowtie` is another option), we need to index the reference genome. This can be done with `bwa index`.

```
bwa index REL606.fa
```

Now, let's map the reads to the reference genome. This is also a fairly intensive step that may take a few minutes.

```
bwa aln REL606.fa SRR098038.fastq.gz > SRR098038.sai
```



![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1 - Look within the ftp directories for these bacterial genome projects. Are there additional genomes (.fna files) present and, if so, what might they represent?


## Annotating both genomes

* Next, annotate both genomes using `prokka`

```
#this may take a while so be patient... Remember, these are full bacterial genomes as opposed to small mitochondrial contigs
prokka K12.fna --outdir K12 --norrna --notrna
prokka O157H7.fna --outdir O157H7 --norrna --notrna
```

## Generating gene lists

* Next, make text files of the predicted gene lists for each genome

```
cd K12/
#generate a gene list text file by grepping the gene names from the .tbl file
cat PROKKA*.tbl | awk '{if ($1 == "gene") {print $2}}' | awk -F'_' '{print $1}' | sort > genelist_K12.txt

cd ../O157H7/
cat PROKKA*.tbl | awk '{if ($1 == "gene") {print $2}}' | awk -F'_' '{print $1}' | sort > genelist_O157H7.txt
```


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2 - How many genes are present in each genome?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3 - Can you find any gene duplicates? If so, provide an example.

## Comparing gene lists

Now, let's compare these lists to find genes that are common to both, duplicated or unique to either

A lot of this work can be done simply using the `comm` and `uniq` functions within the shell.
Explore the command-line usage and options of these commands using `man`

### Comparison including gene duplicates

We can compare both gene lists like this:

```
cd ../  #go back one folder
comm O157H7/genelist_O157H7.txt K12/genelist_K12.txt >geneListComparison.txt
```

Examine the output of `geneListComparison.txt` using `less`

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q4 - What do the genes in column 1, column 2, and column 3 represent? 

Now, suppose we want to output the genes in column 1 (ignoring spaces). We can do so like this:

```
cat geneListComparison.txt | awk -F'\t' '{print $1}' | grep -v -e '^$'
```

... and we can then count them by piping this command to `wc -l`

```
cat geneListComparison.txt | awk -F'\t' '{print $1}' | grep -v -e '^$' | wc -l
```

* Analyze the core versus variable gene content for these two strains

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q5 - Construct a simple Venn diagram (by hand or drawing tool) illustrating the number of shared vs genome-specific genes


### Comparison without gene duplicates (finding unique genes)

In the Venn diagram you have just constructed, if a gene G is duplicated (2 copies) in genome A and not in genome B (1 copy), then gene G will contribute to the unique genes in genome A. This is not really true, since genome B has a copy of gene G.

Let's account for this by "removing gene redundancy". This can be done by first filtering each initial gene list using the tool `uniq`

```
cd K12/
uniq genelist_K12.txt > unique_genelist_K12.txt

cd ../O157H7/
uniq genelist_O157H7.txt > unique_genelist_O157H7.txt

cd ..
```

Now, when we compare these lists using `comm`, we will only be comparing single copies of each gene. Therefore, the result of `comm` should show us only those genes that are unique to genome 1, unique to genome 2, or shared between both

```
comm O157H7/unique_genelist_O157H7.txt K12/unique_genelist_K12.txt > uniqueGeneListComparison.txt
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q6 - Construct a second Venn diagram illustrating the number of shared vs genome-specific genes


## Going further: inspecting your duplicated and unique genes within each organism

Now, let's examine the output of the lists above.

Examine `geneListComparison.txt` to find the gene expansions specific to enterohemorrhagic E. coli O157:H7. The following code will sort the O157:H7-specific genes by their copy number. This will identify those that have undergone the most pathogen-specific duplication.

```
cat geneListComparison.txt | awk -F'\t' '{print $1}' | sort | uniq -c | sort -
n -r | head -20
```

Examine your result carefully. Column 1 states the copy number and copy 2 states the gene name.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q7 - Find a BIOLOGICALLY INTERESTING result among the top-ranked (most duplicated) genes in your list and explain why YOU THINK it is interesting. There are many. Feel free to use Google to help you.

* Repeat the above analysis but using your `uniqueGeneListComparison.txt` file.

Examine the list of genes that are unique to O157H7 (make sure this list is sorted alphabetically using `sort`).

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q8 - Can you find an operon (this will often be a series of genes with similar names) that is unique to E. coli O157H7 that may play a role in its virulence/pathogenicity? Again, use Google to help you.


---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) above.
Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.









 












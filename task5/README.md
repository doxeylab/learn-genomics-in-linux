# Task5 - Comparative genomics - gene set comparison

This task is a tutorial on comparative genomics with a focus on gene set comparison.

You are going to download, annotate and compare the genomes of a enterohemorrhagic <i>E. coli</i> (strain [O157:H7](https://en.wikipedia.org/wiki/Escherichia_coli_O157:H7)) versus a non-pathogenic <i>E. coli</i> (strain [K12](https://en.wikipedia.org/wiki/Escherichia_coli_in_molecular_biology#K-12)). You are then going to identify genes and gene duplications that are unique to each organism and biologically interpret your results.


### Requirements

* Access to a linux-based OS running BASH
* [BLAST](http://blast.ncbi.nlm.nih.gov/)
* [prokka](https://github.com/tseemann/prokka)


## Getting Started

* Login to your linux environment and create a new folder for your task5.
* Work on your assignment in the folder your created.



## Retrieving the raw data

You will be comparing two genomes of <i>E. coli</i> - strain K12 (non-pathogenic lab strain) and O157H7 (pathogenic <i>E. coli</i> associated with disease outbreaks).

* Download both of these genomes using `wget`. Check the man pages for `wget` or use the --help flag to determine how to save the files with the following file names. 
  * Name the O157H7 genome O157H7.fna
  * Name the K12 genome K12.fna

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Escherichia_coli_O157H7_EDL933_uid259/AE005174.fna
ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid20079/CP000948.fna
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1 - Look within the ftp directories for these bacterial genome projects. What do the other files contain (i.e., go to: https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Escherichia_coli_O157H7_EDL933_uid259)?

## Annotating both genomes

* Next, annotate both genomes using `prokka`.

```
#this may take a while so be patient... Remember, these are full bacterial genomes as opposed to small mitochondrial contigs.
prokka O157H7.fna --outdir O157H7 --norrna --notrna
prokka K12.fna --outdir K12 --norrna --notrna
```

## Generating gene lists

* Next, make text files of the predicted gene lists for each genome. Have a look back at Task1 for how to redirect program output to a file.
    * Redirect the output for O157H7 to a file named genelist_O157H7.txt inside the O157H7 folder.
    * Redirect the output for K12 to a file named genelist_K12.txt inside the K12 folder.

```
#generate a gene list text file by grepping the gene names from the .tbl file
cat PROKKA*.tbl | awk '{if ($1 == "gene") {print $2}}' | awk -F '_' '{print $1}' | sort 
```


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2 - How many genes are present in each genome?


## Comparing gene lists

Now, let's compare these lists to find genes that are common to both, duplicated or unique to either.

A lot of this work can be done simply using the `comm` and `uniq` functions within the shell.
Explore the command-line usage and options of these commands using `man`.

### Comparison including gene duplicates

We can compare both gene lists like this in your task5 folder:

```
comm O157H7/genelist_O157H7.txt K12/genelist_K12.txt >geneListComparison.txt
```

Examine the output of `geneListComparison.txt` using `less`.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q3 - What do the genes in column 1, column 2, and column 3 represent? 

Now, suppose we want to output the genes in column 1 (ignoring spaces). We can do so like this:

```
cat geneListComparison.txt | awk -F '\t' '{print $1}' | grep -v -e '^$'
```

... and we can then count them by piping this command to `wc -l`.

```
cat geneListComparison.txt | awk -F '\t' '{print $1}' | grep -v -e '^$' | wc -l
```

* Analyze the core versus variable gene content for these two strains.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q4 - How many genes are only in the O157H7 genome? Only in the K12 genome? In both?


### Comparison without gene duplicates (finding unique genes)
Construct a simple Venn diagram  illustrating the number of shared vs genome-specific genes. In this Venn diagram, if a gene G is duplicated (2 copies) in genome A and not in genome B (1 copy), then gene G will contribute to the unique genes in genome A. This is not really true, since genome B has a copy of gene G.

Let's account for this by "removing gene redundancy". This can be done by first filtering each initial gene list using the tool `uniq`.

```
cd /O157H7
uniq genelist_O157H7.txt > unique_genelist_O157H7.txt

cd ../K12
uniq genelist_K12.txt > unique_genelist_K12.txt

cd ..
```

Now, when we compare these lists using `comm`, we will only be comparing single copies of each gene. Therefore, the result of `comm` should show us only those genes that are unique to genome 1, unique to genome 2, or shared between both

```
comm O157H7/unique_genelist_O157H7.txt K12/unique_genelist_K12.txt > uniqueGeneListComparison.txt
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q5 - How many unique genes are only in the O157H7 genome? Only in the K12 genome? In both?


## Going further: inspecting your duplicated and unique genes within each organism

Now, let's examine the output of the lists above.

Examine `geneListComparison.txt` to find the gene expansions specific to enterohemorrhagic <i>E. coli</i> O157:H7. The following code will sort the O157:H7-specific genes by their copy number. This will identify those that have undergone the most pathogen-specific duplication.

```
cat geneListComparison.txt | awk -F '\t' '{print $1}' | sort | uniq -c | sort -n -r | head -20
```

Examine your result carefully. Column 1 states the copy number and column 2 states the gene name.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q6 - Which gene in O157H7 occurs most frequently? Which genes in K12 occur most frequently?



---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) above. Please submit your answers using the quiz on LEARN.

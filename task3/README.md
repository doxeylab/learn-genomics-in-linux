# Task3 - Genome Annotation

This task is a tutorial on genome annotation using `prokka` and other tools.

You will learn how to perform basic genome annotation, and also how to extract specific regions of interest from your genome sequence.

### Requirements

* Access to a linux-based OS running BASH
* [BLAST](http://blast.ncbi.nlm.nih.gov/)
* [Prokka](https://github.com/tseemann/prokka)
* [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/)
* [eggnog-mapper](http://eggnogdb.embl.de/#/app/emapper)


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

## Annotating your genome from Task2 using `prokka`

By marking the ORFs in your genome (given a min size threshold), you have essentially performed a simple gene finding algorithm. However, there are more advanced ways of gene-finding that take additional criteria into account.

A popular genome annotation tool for prokaryotic genomes is [`prokka`](https://github.com/tseemann/prokka).
`prokka` automates a series of genome annotation tools and is simple to run. It has been installed for you on the server.

* Run prokka using the following command

```
prokka abyss-assembly-contigs.fa
```
<b>Note: This will generate a folder called PROKKA-XXXXXXXX where XXXXXXXX is the current date. It will be different for you than in the examples below.</b>

* Now, locate and download the .gbk file that was produced and view it in `artemis`

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1) Why are there vertical black lines in the middle of predicted ORFs?

To correct this issue, re-start artemis and change your artemis 'Options' to better reflect the source of this genome.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2) What source organism did you select to correct the gene annotations? Paste a screenshot of the artemis genome visualization.

When `prokka` is run without any parameters, it selects 'bacteria' as the default taxonomy.

Look at the `--kingdom` options in `prokka -h` and re-run prokka to use the correct annotation mode.

* Again, open your .gbk file in artemis

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3) Has anything changed in this genome annotation? Provide an example. 


## Annotation of an <i>E. coli</i> genome using `prokka`

Next, let's perform genome annotation on a larger scale.

* Download (or copy from your task1 folder) the E. coli K12 genome from task1 and run `prokka`

```
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/e-coli-k12-genome.fasta.gz
gunzip e-coli-k12-genome.fasta.gz
prokka e-coli-k12-genome.fasta
```

Next, explore the files produced by `prokka`. Start with the .txt file.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q4) How many genes, rRNAs, tRNAs, and CRISPR loci were predicted? What is the size of the genome in Mb?

Prokka also annotates genes based on [COGs](https://www.ncbi.nlm.nih.gov/COG/) and also [E.C.](https://enzyme.expasy.org/) (enzyme commission) numbers. This information can be found in the .tbl file. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q5) How many genes were annotated with COGS? What proportion of total genes is this?

Column 6 of this .tsv file lists the COGs. To print out only column 6, you can use the `cut` command as follows:

```
cut -f6 PROKKA_09182018.tsv
```

Next, let's sort this list using `sort` and pipe it to the `uniq` command.

```
cut -f6 PROKKA_09182018.tsv | sort | uniq
```

How many unique COGs were assigned? Count the number of lines in the output like this:

```
cut -f6 PROKKA_09182018.tsv | sort | uniq | wc -l
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q6) How many unique enzymatic activities (E.C. numbers) were assigned to the E. coli genome?


Next, download the .gbk file produced by prokka to your local machine and view it in `artemis`.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q7) View the genome in DnaPlotter (under File menu) and paste a screenshot.



## Assigning GO terms

Next, we will be assigning Gene Ontology ([GO](http://geneontology.org/)) terms to your predicted genes/proteins.

Prokka identifies homologs of your proteins within the UniProtKB database. Since there are already pre-computed GO terms for all proteins in UniProtKB, we can map these GO terms over using the following commands:

```

```

This will generate an `.GOannotations` file, which contains your predicted functional annotations.

This one-liner will extract column 6 (GO terms), and list them according to their frequency in your proteome.

```
cat .annotations | grep -v "#" | cut -f6 | tr , '\n' | sort | uniq -c | sort -n -r
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q8) What is the most common GO term (GO ID and its function) and why do you think this term is so common? Note: you can get the GO description here: http://amigo.geneontology.org/amigo/term/GO:XXXXXXXX


## After annotation: Extracting genes and regions of interest

Once your genome has been assembled and annotated, you may be interested in identifying and extracting specific genes or regions of interest.

### Extracting genes of interest
For example, suppose you are interested in the "trpE" gene from E. coli. You can see whether this gene exists in the predictions like this:

```
grep "trpE" PROKKA_09202018.tsv
```

This will output
>JKMANJED_01263  CDS     1563    trpE    4.1.3.27        COG0147 Anthranilate synthase component 1

... which tells you that "trpE" has been assigned to the gene labeled "JKMANJED_01263"

You can then extract this gene seqeunce from the gene predictions file (.ffn) like this:

```
# index the .ffn file so we can extract from it
makeblastdb -in PROKKA_09202018.ffn -dbtype 'nucl' -parse_seqids

blastdbcmd -entry JKMANJED_01263 -db PROKKA_09202018.ffn
```

This will output:

>\>>JKMANJED_01263 Anthranilate synthase component 1
ATGCAAACACAAAAACCGACTCTCGAACTGCTAACCTGCGAAGGCGCTTATCGCGACAATCCCACCGCGCTTTTTCACCA
GTTGTGTGGGGATCGTCCGGCAACGCTGCTGCTGGAATCCGCAGATATCGACAGCAAAGATGATTTAAAAAGCCTGCTGC
TGGTAGACAGTGCGCTGCGCATTACAGCTTTAGGTGACACTGTCACAATCCAGGCACTTTCCGGCAACGGCGAAGCCCTC
CTGGCACTACTGGATAACGCCCTGCCTGCGGGTGTGGAAAGTGAACAATCACCAAACTGCCGTGTGCTGCGCTTCCCCCC
...

<b>Note: in this example we searched the annotations with a text query "trpE". However, the best way of finding your gene of interest is to do a BLAST search since it may not be labeled correctly in your annotations</b>

### Extracting regions of interest

Next, suppose are interested in extracting the promoter of the "trp operon". The trp operon promoter can be found directly upstream of the <b>trpE</b> gene. These regions are not in the annotations files so you will need to locate them yourself.

First, let's see where the trpE gene is located in the genome:

```
grep "trpE" PROKKA_09182018.gff
```

This will output:

>U00096.3        Prodigal:2.6    CDS     1321384 1322946 .       -       0       ID=JKMANJED_01263;eC_number=4.1.3.27;Name=trpE;db_xref=COG:COG0147;gene=trpE;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P00895;locus_tag=JKMANJED_01263;product=Anthranilate synthase component 1

This tells us that trpE is located in entry "U00096.3" at chromosome position "1321384 to 1322946" and encoded on the minus (-) strand.

To extract the sequence for these coordinates, we can use `blastdbcmd` against the genome as follows:

```
# index the genome so we can extract regions from it
makeblastdb -in PROKKA_09182018.fna -dbtype 'nucl' -parse_seqids

blastdbcmd -entry U00096.3 -db PROKKA_09182018.fna -range 1321384-1322946 -strand minus
```

This should produce a FASTA sequence output of the gene identical to that in the above example.

But you are not interested in the gene sequence; you actually want the promoter region.

* So, modify the code above to extract the 30-nucleotide long sequence upstream of the trpE gene

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q8) Paste this sequence into your assignment as well as the code you used to extract it.


### Advanced: Extracting the rRNAs predicted by barrnap

Sometimes you may be interested in extracting multiple genes or regions at once. E.g., suppose you want to extract all of the regions corresponding to predicted 16S rRNA sequences. In `prokka`, rRNA genes are predicted for you using the `barrnap` tool.

Here is a two-liner to extract the 16S rRNAs predicted by `barrnap`, and for fun we will also pipe this to `muscle` to do an automatic alignment.

```
cat *.gff | grep "barrnap" | awk '{ if ($7 == "-") {print $1" "$4"-"$5" minus"} else {print $1" "$4"-"$5" plus"} }' >rRNAs.txt
blastdbcmd -db PROKKA_09182018.fna -entry_batch rRNAs.txt | muscle
```


## Analyzing a mystery genome of unknown source

And now for something a little more difficult.

* Download this "mystery" genome of unknown source

https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task3/mysteryGenome.fna.gz

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q9) Identify a full-length 16S rRNA sequence. Paste this sequence into your assignment and include your source code.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q10) Now, using [web-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn), search this sequence against the NCBI 16S database. What is the taxonomic origin of this genome? 

![](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/task3/16Ssearch.png)

---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) above.
Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.



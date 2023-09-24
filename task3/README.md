# Task3 - Genome Annotation

This task is a tutorial on genome annotation using `prokka` and other tools.

You will learn how to perform basic genome annotation, and also how to extract specific regions of interest from your genome sequence.

### Requirements

Graphical software indicated by (*)

* Access to a linux-based OS running BASH
* [BLAST](http://blast.ncbi.nlm.nih.gov/)
* [Prokka](https://github.com/tseemann/prokka)
* [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) *
Note: for Artemis, you may need to first install the JRE (Java Runtime Environment) and/or JDK (Java Development Kit) on your system. Also, if you have difficulties installing the latest version on Windows, try installing version 17.0.1 here: ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/v17.0.1/.
* uniprot2go.py script located [here](https://github.com/doxeylab/learn-genomics-in-linux/blob/master/task3/uniprot2go.py)
-- must be installed in /usr/bin

## Installation

If you do not already have access to a GUI running the graphical software listed above, please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

---

## Getting Started

* Login to your linux environment as you did in task1.

* Create a new folder for task3.

```
mkdir task3  #creates folder
cd task3 #enters into folder
```

## Retrieving the raw data

* Copy the genome you assembled with `abyss` from task2.

```
cp ../task2/abyss-assembly-contigs.fa . 
```

## Exploring the assembly using `artemis`


Now that we have generated a good quality assembly, let's explore the genome sequence itself and do some very basic annotation using `artemis`. 

Visualize the genome you have produced (using `abyss`) with the `artemis` application. Note: you will need to have `artemis` installed on your local machine.
You will also need to download your contigs to your local machine. Open your contigs.fa file in `artemis`.

- What are the black vertical lines that appear in the sequence window?
- How do you product a gc plot of the genome?
- Why would a researcher create a gc plot?
- How do you mark open reading frames (ORFs)?


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1) How many ORFs are there of length >= 100 amino acids? 

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2) How many ORFs are there of length >= 50 amino acids? 

## Annotating your genome from Task2 using `prokka`

By marking the ORFs in your genome (given a min size threshold), you have essentially performed a simple gene finding algorithm. However, there are more advanced ways of gene-finding that take additional criteria into account.

A popular genome annotation tool for prokaryotic genomes is [`prokka`](https://github.com/tseemann/prokka).
`prokka` automates a series of genome annotation tools and is simple to run. It has been installed for you on the server.

* Run `prokka` using the following command.

```
prokka abyss-assembly-contigs.fa
```
<b>Note: This will generate a folder called PROKKA-XXXXXXXX where XXXXXXXX is the current date. It will be different for you than in the examples below.</b>

* Now, locate and download the .gbk file that was produced and view it in `artemis`.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q3) You will notice that there are vertical black lines in the middle of predicted ORFs. What do these lines represent?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q4) Re-start `artemis` and change your `artemis` 'Options' to better reflect the source of this genome. Which source did you choose (e.g., "Standard", "Vertebrate Mitochondrial", etc.)?

When `prokka` is run without any parameters, it selects 'bacteria' as the default taxonomy.

Look at the `--kingdom` options in `prokka -h` and re-run `prokka` to use the correct annotation mode. You will also need to use `--outdir` to specify a folder for your new results.

* Again, open your .gbk file in `artemis`.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q5) Has anything changed in this genome annotation? Examine the CDSs, tRNAs, and rRNAs, and their annotations.


## Annotation of an <i>E. coli</i> genome using `prokka`

Next, let's perform genome annotation on a larger scale.

* Download (or copy from your task1 folder) the <i>E. coli</i> K12 genome below from task1 and annotate it using `prokka`

```
https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task1/e-coli-k12-genome.fasta.gz
```

Next, explore the files produced by `prokka`. Start with the .txt file.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q6) How many genes, rRNAs, tRNAs, and CRISPR loci were predicted? What is the size of the genome in Mb?

`prokka` also annotates genes based on [COGs](https://www.ncbi.nlm.nih.gov/COG/) and also [E.C.](https://enzyme.expasy.org/) (enzyme commission) numbers. This information can be found in the .tbl file. 

Column 6 of this .tsv file lists the COGs. To print out only column 6, you can use the `cut` command as follows (replace "yourPROKKAoutput"):

```
cut -f6 yourPROKKAoutput.tsv
```

Using commands such as `cut`, `sort`, `grep`, `uniq`, and `wc` answer the following two questions (Q7 and Q8).

e.g., this line below will count the number of unique entries in column 3 of file.txt

```
cut -f3 file.txt | sort | uniq | wc -l
```


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q7) How many genes were annotated with COGs?


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q8) How many unique enzymatic activities (E.C. numbers) were assigned to the <i>E. coli</i> genome? Note: `1.-.-.-` and `1.1.1.17` would count as two separate E.C. numbers.



## Assigning GO terms

Next, we will be assigning Gene Ontology ([GO](http://geneontology.org/)) terms to your predicted genes/proteins.

`prokka` identifies homologs of your proteins within the UniProtKB database. Since there are already pre-computed GO terms for all proteins in UniProtKB, we can map these GO terms over using the following commands:

```
#extract the predicted proteins that have been mapped to entries in UniProt
cat yourPROKKAoutput.gff | grep -o "UniProtKB.*;" | awk -F'[:;=]' '{print $4" "$2}' >uniProts.txt

#assign GO annotations from a uniprot-GO database table
python2.7 /usr/bin/uniprot2go.py -i uniProts.txt -d /fsys1/data/uniprot2go/uniprot-vs-go-db.sl3 >go.annotations
```

This will generate an `go.annotations` file, which contains your predicted functional annotations.

This one-liner will extract column 3 (GO terms), and list the top 20 according to their frequency in your proteome.

```
cat go.annotations | awk '{print $3}' | tr "," "\n" | sort | uniq -c | sort -n -r | head -20
```

Now, there is a lot you can explore using your predicted GO terms for your genome.
e.g., Suppose you want to find all the predicted DNA binding proteins. Look [here](http://amigo.geneontology.org/amigo) to find the GO ID for "DNA binding".

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q9) How many proteins were annotated with the GO term for "DNA binding"? 

## After annotation: Extracting genes and regions of interest

Once your genome has been assembled and annotated, you may be interested in identifying and extracting specific genes or regions of interest.

### Extracting genes of interest
For example, suppose you are interested in the "trpE" gene from <i>E. coli</i>. You can see whether this gene exists in the predictions like this:

```
grep "trpE" yourPROKKAoutput.tsv
```

This will output
>JKMANJED_01263  CDS     1563    trpE    4.1.3.27        COG0147 Anthranilate synthase component 1

... which tells you that "trpE" has been assigned to the gene labeled "JKMANJED_01263".

You can then extract this gene sequence from the gene predictions file (.ffn) like this:

```
# index the .ffn file so we can extract from it
makeblastdb -in yourPROKKAoutput.ffn -dbtype 'nucl' -parse_seqids

blastdbcmd -entry JKMANJED_01263 -db yourPROKKAoutput.ffn
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

Next, suppose you are interested in extracting the promoter of the "trp operon". See [here](https://en.wikipedia.org/wiki/Trp_operon) for some background information. The trp operon regulatory sequences (operator and promoter) can be found upstream of the <b>trpE</b> gene. These regions are not in the annotations files so you will need to locate them yourself.

First, let's see where the trpE gene is located in the genome:

```
grep "trpE" yourPROKKAoutput.gff
```

This will output:

>U00096.3        Prodigal:2.6    CDS     1321384 1322946 .       -       0       ID=JKMANJED_01263;eC_number=4.1.3.27;Name=trpE;db_xref=COG:COG0147;gene=trpE;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P00895;locus_tag=JKMANJED_01263;product=Anthranilate synthase component 1

This tells us that trpE is located in entry "U00096.3" at chromosome position "1321384 to 1322946" and encoded on the minus (-) strand.

To extract the sequence for these coordinates, we can use `blastdbcmd` against the genome as follows:

```
# index the genome so we can extract regions from it
makeblastdb -in yourPROKKAoutput.fna -dbtype 'nucl' -parse_seqids

blastdbcmd -entry U00096.3 -db yourPROKKAoutput.fna -range 1321384-1322946 -strand minus
```

This should produce a FASTA sequence output of the gene identical to that in the above example.

But you are not interested in the gene sequence; you actually want the upstream regulatory region. Suppose you want to identify the 30-nucleotide long region upstream (before but not including the start codon) of the trpE coding sequence. By modifying the code above, answer the following question.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q10) What is the 30-nucleotide long sequence immediately upstream of the TrpE coding sequence?


### Extracting the rRNAs predicted by barrnap

Sometimes you may be interested in extracting multiple genes or regions at once. E.g., suppose you want to extract all of the regions corresponding to predicted 16S rRNA sequences. In `prokka`, rRNA genes are predicted for you using the `barrnap` tool.

Here is a two-liner to extract the 16S rRNAs predicted by `barrnap`.

```
cat yourPROKKAoutput.gff | grep "barrnap" | awk '{ if ($7 == "-") {print $1" "$4"-"$5" minus"} else {print $1" "$4"-"$5" plus"} }' >rRNAs.txt
blastdbcmd -db yourPROKKAoutput.fna -entry_batch rRNAs.txt > rRNAs.fa
```

Now, to predict taxonomy, we can BLAST these rRNA sequences against the NCBI 16S database for example using [web-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn). Note, that there may be multiple rRNAs and some of them may be partial sequences.

![](https://github.com/doxeylab/learn-genomics-in-linux/blob/master/task3/16Ssearch.png)


## Analyzing a mystery genome of unknown source

And now for something a little more difficult.

* Download this "mystery" genome of unknown source.

https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task3/mysteryGenome.fna.gz

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q11) Based on 16S rRNA sequences, what is the taxonomic origin of this genome (genus and species)? 
e.g., "<i>Escherichia coli</i>"


---

# ASSIGNMENT QUESTIONS


Please answer questions 1-11 above on LEARN under Quizzes.


#

Congratulations. You have now completed Task 3.



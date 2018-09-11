# Task2 - Genome Assembly

In this task, you will download raw sequencing data, perform genome assembly, genome annotation, and then compare your assembly to the database using BLAST.

### Requirements

* Access to a linux-based OS running BASH
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
* [abyss](https://github.com/bcgsc/abyss)
* [tablet](https://ics.hutton.ac.uk/tablet/)
* [artemis](http://sanger-pathogens.github.io/Artemis/Artemis/)

---

## Getting Started

Login to your linux environment as you did in task1. e.g.,

```
ssh -i /path/to/your/.ssh/publickey yourUserName@remoteIP
```

Create a new folder for your task2

```
mkdir task2  #creates folder
cd task2 #enters into folder
```

## Retrieving the raw data

There are several ways to download data. Two common tools are `curl` and `wget`.

The following command will download the raw sequencing data into your folder and then uncompress it.

```
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task2/mt_reads.fastq.gz

gunzip mt_reads.fastq.gz

```

Now explore the structure of the fastq file

```
less mt_reads.fastq.gz
```

:large_blue_circle: <b>Q1</b> Paste the lines corresponding to the first 3 reads.


:large_blue_circle: <b>Q2</b> How many lines are in the file? Hint: you can use `grep` or `wc` to answer this


## Data preprocessing

Before we can assemble a genome, we need to:

* Assess the quality of the sequencing data
* Demultiplex the data
* Trim barcodes
* Filter out low-quality reads (this is called quality filtering)

### Quality assessment
For a quick quality report, you can use the program `fastqc`.
The command below will analyze the mt_reads.fastq file and place the results in the 'qcresults' folder.

```
fastqc -o qcresults mt_reads.fastq
```

Now open the HTML file that was generated in firefox. The ampersand (&) will run the command in the background so you can still use the terminal.

```
firefox qcresults/mt_reads_fastqc/fastqc_report.html &
```

### Splitting the barcodes (demultiplexing)

Sequencing data may be barcoded. In this case, the data contains two different samples, each with a unique barcode.
This allows us to split the data by sample. Sometimes, sequencing data can have tens or hundreds of barcodes. See [multiplexing](https://www.illumina.com/science/technology/next-generation-sequencing/multiplex-sequencing.html)

We will use a standard script from the `fastx toolkit` to split the data by its known barcodes (defined already for you in the file downloaded below)

```
#first download the barcodes file
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task2/mt_barcodes.txt

#now split
fastx_barcode_splitter.pl <mt_reads.fastq --bcfile mt_barcodes.txt --bol --suffix .fastq --prefix results/
```

There are now two .fastq files in the results directory; one for each barcode.  There is also an unmatched.fasta file which should be empty.  We will be focusing on the first sample, ie. the one now in results/mt1.fastq.

### Barcode trimming

Barcode trimming is needed to remove the barcode sequences from the beginning of each read. The Q33 is required due to differences in sanger and illumina encoding.

```
cd results
fastx_trimmer -i mt1.fastq -f 8 -o trimmed_mt1.fastq -Q33
```


### Quality filtering

Next, we need to remove low quality sequences. This will increase the accuracy of the assembly.

```
fastq_quality_filter -i trimmed_mt1.fastq -q 25 -p 80 -o qual_trim_mt1.fastq -Q33 -v
```


## Genome Assembly

Now we are ready to assemble a genome. 

### Assembly with Velvet
To start we are going to try using the popular `velvet` assembler: https://en.wikipedia.org/wiki/Velvet_assembler
Like many assemblers, `velvet` performs genome assembly using de bruijn graphs. This means that we must choose a value of <b>k</b> to define the k-mers (sequence fragments of length k) to be used in constructing the graph.

https://en.wikipedia.org/wiki/De_Bruijn_graph
https://en.wikipedia.org/wiki/De_novo_sequence_assemblers



### Assembly with ABYSS

Inspect your contigs.fa file (should be in your velvet_k12 directory). There should be 139 separate contigs within this file. You can verify this by running:

```
grep -c ">" contigs.fa
```

But 139 separate contigs is not ideal. To improve the assembly, try assembling with the popular abyss assembler (with a kmer length of 21) using the following command:

```
abyss-pe k=21 in='qual_trim_mt1.fastq' name=abyss-contigs
```

:large_blue_circle: How many contigs are produced now?


Examine the dinucleotide content of your assembled mitochondrial genome using the tool compseq. 

:large_blue_circle: What dinucleotide is most underrepresented in the genome? :grey_question:


## Visualization using Artemis

```
art abyss_contigs.fa
```

:large_blue_circle: Produce a gc plot (& increase sliding window length to maximum) covering the whole mitochondrial genome and paste a screenshot.


## Genome Annotation

...

## BLAST



# ASSIGNMENT QUESTIONS

The questions for Task 2 are indicated by the lines starting with :large_blue_circle: above
Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.

# Task2 - Genome Assembly

In this task, you will download raw sequencing data, perform genome assembly, genome annotation, and then compare your assembly to the database using BLAST.

### Requirements

* Access to a linux-based OS running BASH
* <b>This task also requires graphical software indicated below </b> (*)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) *
* [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
* [abyss](https://github.com/bcgsc/abyss)
* [tablet](https://ics.hutton.ac.uk/tablet/) *
* [artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) *
* [bandage](http://rrwick.github.io/Bandage/) *

## Installation

If you do not already have access to a GUI running the graphical software listed above, please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

All software used are available for Mac/Windows/Linux.

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

Download the raw sequencing data into your folder and then uncompress it.

```
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task2/mt_reads.fastq.gz

gunzip mt_reads.fastq.gz

```

Now explore the structure of the fastq file

```
less mt_reads.fastq.gz
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1) Paste the lines corresponding to the first 3 reads.


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2) How many reads are in the file? Hint: you can use `grep` to answer this.


## Data preprocessing

Before we can assemble a genome, we need to:

* Assess the quality of the sequencing data
* Demultiplex the data
* Trim barcodes
* Filter out low-quality reads (this is called quality filtering)

### Quality assessment
For a quick quality report, you can use the program `fastqc`.
The command below will analyze the mt_reads.fastq file and produce an .html results file.

```
fastqc mt_reads.fastq
```

Next, [<b>download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) the 'fastqc_report.html' file to your local machine and open it in a web browser.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3) Inspect the quality of your reads by examining the k-mer abundance profile. Why do you think some k-mers are unusually abundant?


### Splitting the barcodes (demultiplexing)

Sequencing data may be barcoded. In this case, the data contains two different samples, each with a unique barcode.
This allows us to split the data by sample. Sometimes, sequencing data can have tens or hundreds of barcodes. See [multiplexing](https://www.illumina.com/science/technology/next-generation-sequencing/multiplex-sequencing.html)

We will use a standard script from the `fastx toolkit` to split the data by its known barcodes (defined already for you in the file downloaded below)

```
#first download the barcodes file
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task2/mt_barcodes.txt

#now split
fastx_barcode_splitter.pl <mt_reads.fastq --bcfile mt_barcodes.txt --bol --suffix .fastq --prefix splitData_
```

There are now two .fastq files; one for each barcode.  There is also an unmatched.fasta file which should be empty.  We will be focusing on the first sample, ie. the one now in the file 'splitData_mt1.fastq'.

### Barcode trimming

Barcode trimming is needed to remove the barcode sequences from the beginning of each read. The Q33 is required due to differences in sanger and illumina encoding.

```
fastx_trimmer -i splitData_mt1.fastq -f 8 -o trimmed_mt1.fastq -Q33
```


### Quality filtering

Next, we need to remove low quality sequences. This will increase the accuracy of the assembly.

```
fastq_quality_filter -i trimmed_mt1.fastq -q 25 -p 80 -o qual_trim_mt1.fastq -Q33 -v
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q4) Examine the quality of your data again using `fastqc`. Has the quality improved? Justify your answer.


## Genome Assembly

Now we are ready to assemble a genome. 

### Assembly with Velvet
To start we are going to try using the popular `velvet` assembler: https://en.wikipedia.org/wiki/Velvet_assembler
Like many assemblers, `velvet` performs genome assembly using de bruijn graphs. This means that we must choose a value of <b>k</b> to define the k-mers (sequence fragments of length k) to be used in constructing the graph.

https://en.wikipedia.org/wiki/De_Bruijn_graph
https://en.wikipedia.org/wiki/De_novo_sequence_assemblers


The command below will compute the graph. The first parameter is the folder name (you choose) and the second parameter is the value of k.

```
velveth out_21 21 -short -fastq trimmed_mt1.fastq
```

Now, to compute the actual contig sequences from the graph, run the following:

```
velvetg out_21/ -scaffolding no -read_trkg yes -amos_file yes
```

Inspect the contigs.fa file that has been produced (will be in out_21 folder).

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q5) How many contigs did you produce? Try varying the k-mer value. Did the number of contigs change?

### Assembly visualization with Tablet

Velvet has the option of keeping track of where the reads map to the assembly using the `-read_trkg` flag. This will produce a `velvet_asm.afg` file.

[<b>Download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) this file to your local machine. Tip: find the path to your file with `realpath yourFile.txt`

Then open it in `tablet`.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q6) Report the average contig length and N50
![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q7)  Choose any contig and view its assembly using tablet. Paste a screenshot.

### Assembly visualization with Bandage

Velvet and other de bruijn assemblers produce a graph that can be visualized. `bandage` is an excellent tool for this purpose.

Find the 'lastgraph' file produced by `velvet` and [<b>download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) it to your local machine. Tip: find the path to your file with `realpath yourFile.txt`

Open this file in the `Bandage` application. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q8) How long is the largest connected component of this graph?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q9) Copy and paste a screenshot of a region of the graph that has an ambiguous assembly.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q10) Copy and paste a screenshot of a region that has an unambiguous assembly (region should contain more than one node).

### Assembly with ABYSS

Next, let's try with a different assembler. We will be using the popular [abyss](https://github.com/bcgsc/abyss) assembler. We will keep the value of k = 21.

```
abyss-pe k=21 in='qual_trim_mt1.fastq' name=abyss-assembly
```

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q11) How many contigs did abyss generate?

Take a look at your abyss-assembly-stats output file.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q12) How long (in kB) is your assembly?

#### What is the taxonomic source of your genome? Explore with BLAST

You still do not know the source of this genome. Is it eukaryotic? bacterial? is it nuclear or mitochondrial?

To investigate this question, do a BLAST search using the <b>online</b> [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) tool.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q13) Based on the BLAST result, what organism do you think this genome came from and what kind of genome is it?

## Genome visualization using Artemis

Next, let's visualize the contigs you have produced using Artemis. Note: you will need to have `artemis` installed on your local machine.
You will also need to [<b>downloaded</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) your contigs to your local machine. Tip: find the path to your file with `realpath yourFile.txt`

Open your contigs.fa file in `artemis`.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q14) What are the black vertical lines that appear in the sequence window?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q15) Next, produce a gc plot of the genome (look under the 'Graph' window). Increase the sliding window length to 500 and make sure the whole genome is visible. Paste a screenshot.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q16) Now, mark the open reading frames (ORFs) with min length = 100nt (under 'Create' menu). Paste a snapshot of your Artemis window (make sure the full genome is visible). How many ORFs of length > 100nt are there? How many ORFs of length > 50nt are there? 


## Genome Annotation

By marking the ORFs in your genome (given a min size threshold), you have essentially performed a simple gene finding algorithm. However, there are more advanced ways of gene-finding that take additional criteria into account.

A popular genome annotation tool for prokaryotic genomes is [`prokka`](https://github.com/tseemann/prokka).
`prokka` automates a series of genome annotation tools and is simple to run. It has been installed for you on the server.

Type

```
prokka abyss-assembly-contigs.fa
```

* Now, download the .gbk file that was produced and view it in Artemis

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q17) Are the predicted gene locations consistent with your earlier ORF predictions? 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q18) Why are there vertical black lines in the middle of predicted ORFs?

* Quit artemis and change your artemis 'Options' to better reflect the source of this genome. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q19) Does this correct the issues above? Paste a screenshot.



# ASSIGNMENT QUESTIONS

The questions for Task 2 are indicated by the lines starting with above
Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.

# Task 2 - Genome Assembly

In this lab, you will download raw sequencing data, perform genome assembly, visualize and analyze your assemblies, and compare the assembled genome sequence to the database using BLAST.

### Requirements

* Access to a linux-based OS running BASH 
* <b>This task also requires graphical software indicated below </b> (*)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
* [abyss](https://github.com/bcgsc/abyss)
* [tablet](https://ics.hutton.ac.uk/tablet/) *
* [bandage](http://rrwick.github.io/Bandage/) *


## Installation

If you do not already have access to a GUI running the graphical software listed above, please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

All software used are available for Mac/Windows/Linux.

---

## Getting Started

Login to your Linux environment as you did in task1, and create a new folder for your task2 work

```
mkdir task2  #creates folder
cd task2 #enters into folder
```

## Retrieving the raw data

Download the raw sequencing data from *https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task2/mt_reads.fastq.gz* into your folder and then uncompress it.

Explore the file using `less`.

Next, consider how would you return the sequence of the nth read using only head and tail or grep? Note that each read consists of a fixed number of lines and the first line of a read starts with a specific character.


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1) What is the sequence of the third read in the file? Make sure to remove all spaces.


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2) How many reads are in the file?  Hint: the command `grep "^X"` reports all lines starting with the character `X`.



## Data preprocessing

Before we can assemble a genome, we need to:

1) Assess the quality of the sequencing data
2) Demultiplex the data
3) Trim barcodes
4) Filter out low-quality reads (this is called quality filtering)

### Quality assessment
For a quick quality report, you can use the program `fastqc`.
The command below will analyze the mt_reads.fastq file and produce an .html results file.

```
fastqc mt_reads.fastq
```

[<b>Download</b>](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task1/gcloud-download.png) the 'fastqc_report.html' file to your local machine and open it in a web browser. Tip: find the path to your file with `realpath yourFile.txt`

Explore and inspect the FastQC report for mt_reads.fastq.


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q3) Which of the following statements is correct?
* The reads passed all of the quality control measures
* The reads failed all of the quality control measures
* The reads passed some of the quality control measures but failed others


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q4) The per-base sequence quality is lowest at the \_\_\_\_\_ of the reads.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q5) Most of the reads were assigned a quality (Phred) score of \_\_\_\_\_ .

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q6) Examine the per-base sequence content. The base composition is unusual/unexpected for position \_\_\_\_\_ to position \_\_\_\_\_ of the reads.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q7) This unexpected composition may be due to the inclusion of \_\_\_\_\_ .




### Splitting the barcodes (demultiplexing)

Sequencing data may be barcoded. In this case, the data contains two different samples, each with a unique barcode.
This allows us to split the data by sample. Sometimes, sequencing data can have tens or hundreds of barcodes. See [multiplexing](https://www.illumina.com/science/technology/next-generation-sequencing/multiplex-sequencing.html)

We will use a standard script from the `fastx toolkit` to split the data by its known barcodes (defined already for you in the file downloaded below)

```
#first download the barcodes file
wget https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task2/mt_barcodes.txt

#now split
fastx_barcode_splitter.pl <mt_reads.fastq --bcfile mt_barcodes.txt --bol --suffix .fastq --prefix splitData_
```

There are now two .fastq files; one for each barcode.  There is also an unmatched.fasta file which should be empty.  We will be focusing on the first sample, ie. the one now in the file 'splitData_mt1.fastq'.



### Barcode trimming

Barcode trimming is needed to remove the barcode sequences from the beginning of each read. The Q33 is required due to differences in sanger and illumina encoding.

```
fastx_trimmer -i splitData_mt1.fastq -f 9 -o trimmed_mt1.fastq -Q33
```


### Quality filtering

Next, we need to remove low quality sequences. This will increase the accuracy of the assembly.

```
fastq_quality_filter -i trimmed_mt1.fastq -q 25 -p 80 -o qual_trim_mt1.fastq -Q33 -v
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q8) What percentage of reads were removed?


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q9) Compare FastQC reports from before and after trimming and quality filtering. Which of the measures improved from a warning/fail to a pass?



## Genome Assembly

Now we are ready to assemble a genome. 

### Assembly with Velvet
To start we are going to try using the popular `velvet` assembler. Like many assemblers, `velvet` performs genome assembly using de bruijn graphs. This means that we must choose a value of <b>k</b> to define the k-mers (sequence fragments of length k) to be used in constructing the graph.

Read more:
[velvet](https://en.wikipedia.org/wiki/Velvet_assembler).
[de bruijn graphs](https://en.wikipedia.org/wiki/De_Bruijn_graph)
[de novo assemblers](https://en.wikipedia.org/wiki/De_novo_sequence_assemblers)


The commands below will compute the graph. The first parameter is the folder name (you choose) and the second parameter is the value of k. So below, we are assembling the genome from the trimmed and quality-filtered reads using a k-mer value of 21.

```
velveth out_21 21 -short -fastq qual_trim_mt1.fastq
```

Next, to compute the actual contig sequences from the graph, run the following:

```
velvetg out_21/ -scaffolding no -read_trkg yes -amos_file yes
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q10) How many nodes are there in the graph that was produced?

Inspect the contigs.fa file that has been produced (will be in out_21 folder).


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q11) How many contigs do you get using k=21?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q12) How many contigs do you get using k=31?


## Assembly visualization

### Assembly visualization with Tablet

Velvet has the option of keeping track of where the reads map to the assembly using the `-read_trkg` flag. This will produce a `velvet_asm.afg` file.

[<b>Download</b>](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task1/gcloud-download.png) this file (from your k=31 assembly) to your local machine. 
Then open it in `tablet`. Tablet is a great program to explore how reads map to assemblies and genomes.


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q13) What is the average contig length?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q14) What is the N50 value?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q15) Examine the read coverage across the longest contig. Does the coverage distribution match that shown  [<b>here</b>](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task2/tablet-coverage-plot.png)?


Explore `tablet` more on your own. We will be using it later in the course.

### Assembly visualization with Bandage (optional)

Velvet and other de bruijn assemblers produce a graph that can be visualized. `bandage` is an excellent tool for this purpose.

If you are interested, locate the 'lastgraph' file produced by `velvet` and [<b>download</b>](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task1/gcloud-download.png) it to your local machine.

Open this file in the `Bandage` application, and explore further. 


### Generating an improved assembly with ABYSS

As you can see based on the results from above, `velvet` (with the parameters we chose) did not yield a high quality assembly. It is too fragmented.

Often in genomics it is useful to try numerous parameters and different assemblers. Let's try to assemble this genome again but with a different assembler. This time, we will be using the popular [abyss](https://github.com/bcgsc/abyss) assembler and we'll keep the value of k = 21. And see [here](https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#stats) for info on stats reported by an Abyss assembly.


```
abyss-pe k=21 in='qual_trim_mt1.fastq' name=abyss-assembly
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q16) How long (# bases) is your assembly?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q17) What is the N50 value (# bases)?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q18) How many contigs did abyss generate?



### What is the taxonomic source of your genome? Explore with BLAST

You still do not know the source of this genome. Is it eukaryotic? bacterial? Is it a nuclear genome, a plasmid, or something else?

To investigate this question, do a BLAST search using the <b>online</b> [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) tool. Use the Abyss assembly as your query.


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q19) Based on the BLAST result, describe the most likely source of this DNA sequence?



# ASSIGNMENT QUESTIONS


Please answer questions 1-19 above on LEARN under Quizzes.


#

Congratulations. You have now completed Task 2.

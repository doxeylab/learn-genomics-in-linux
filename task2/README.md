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

Download the raw sequencing data from *https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task2/mt_reads.fastq.gz* into your folder and then uncompress it.

- How many reads are in the file? Hint: the command `grep "^X"` reports all lines starting with the character `X`.

- How would you return the sequence of the nth read using only head and tail or grep? Hint: each read consists of x lines and the first line of read starts with a specific character.




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

[<b>Download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) the 'fastqc_report.html' file to your local machine and open it in a web browser. Tip: find the path to your file with `realpath yourFile.txt`

- What information do you get from the report?



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
fastx_trimmer -i splitData_mt1.fastq -f 9 -o trimmed_mt1.fastq -Q33
```


### Quality filtering

Next, we need to remove low quality sequences. This will increase the accuracy of the assembly.

```
fastq_quality_filter -i trimmed_mt1.fastq -q 25 -p 80 -o qual_trim_mt1.fastq -Q33 -v
```

- What percentage of reads were removed?
- After creating a new FastQC report what measures have changed?




## Genome Assembly

Now we are ready to assemble a genome. 

### Assembly with Velvet
To start we are going to try using the popular `velvet` assembler. Like many assemblers, `velvet` performs genome assembly using de bruijn graphs. This means that we must choose a value of <b>k</b> to define the k-mers (sequence fragments of length k) to be used in constructing the graph.

Read more:
[velvet](https://en.wikipedia.org/wiki/Velvet_assembler).
[de bruijn graphs](https://en.wikipedia.org/wiki/De_Bruijn_graph)
[de novo assemblers](https://en.wikipedia.org/wiki/De_novo_sequence_assemblers)


The command below will compute the graph. The first parameter is the folder name (you choose) and the second parameter is the value of k.

```
velveth out_21 21 -short -fastq trimmed_mt1.fastq
```

Next, to compute the actual contig sequences from the graph, run the following:

```
velvetg out_21/ -scaffolding no -read_trkg yes -amos_file yes
```

Inspect the contigs.fa file that has been produced (will be in out_21 folder).
- How many contigs do you have?
- How many contigs do you get if you increase the k value to 31?



## Assembly visualization

### Assembly visualization with Tablet

Velvet has the option of keeping track of where the reads map to the assembly using the `-read_trkg` flag. This will produce a `velvet_asm.afg` file.

[<b>Download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) this file to your local machine. 
Then open it in `tablet`. Tablet is a great program to explore how reads map to assemblies and genomes.

- Where can you find the average contig length and N50 value? 


Explore `tablet` more on your own. We will be using it later in the course.

### Assembly visualization with Bandage

Velvet and other de bruijn assemblers produce a graph that can be visualized. `bandage` is an excellent tool for this purpose.

Find the 'lastgraph' file produced by `velvet` and [<b>download</b>](https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task1/gcloud-download.png) it to your local machine.

Open this file in the `Bandage` application. 

- How long is the largest connected component of this graph



### Generating an improved assembly with ABYSS

As you can see based on the results from above, `velvet` (with the parameters we chose) did not yield a high quality assembly. It is too fragmented.

Often in genomics it is useful to try numerous parameters and different assemblers. Let's try to assemble this genome again but with a different assembler. This time, we will be using the popular [abyss](https://github.com/bcgsc/abyss) assembler but we'll keep the value of k = 21.


```
abyss-pe k=21 in='qual_trim_mt1.fastq' name=abyss-assembly
```

- How many contigs did abyss generate?
- According to the abyss-assembly-stats file how long (in kB) is your assembly? What is the N50?
- How many contigs were produced?





#### What is the taxonomic source of your genome? Explore with BLAST

You still do not know the source of this genome. Is it eukaryotic? bacterial? is it nuclear or mitochondrial?

To investigate this question, do a BLAST search using the <b>online</b> [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) tool. Use the Abyss assemble as your query.

-  Based on the BLAST result, what organism do you think this genome came from and what kind of genome is it?




# ASSIGNMENT QUESTIONS
#### Use the file mt_reads.fastq to answer the following questions.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1) What is the sequence of the third read in the file? Make sure to remove all spaces.


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2) How many reads are in the file? 


Generate and inspect the FastQC report for mt_reads.fastq


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3) Which of the following statement is correct?
* The reads passed all of the quality control measures
* The reads failed all of the quality control measures
* The reads passed some of the quality control measures but failed others


#### Use the file splitData_mt2.fastq to answer the following questions.

##### After trimming the barcodes and quality filtering.


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q4) What percentage of reads were removed?


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q5) Compare FastQC reports from before and after trimming and quality filtering. Which of the measures improved from a warning/fail to a pass?


![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q6) How many contigs do you get when assembling using Velvet with a k=31?


##### Download the velvet_asm.afg file and open it in tablet.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q7) What is the average contig length?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q8) What is the N50 value?

##### Download the lastgraph file and open with Bandage

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q9) How long is the largest connected component of the graph?

##### Run Abyss with a k=21 and use the abyss-assembly-stats to answer the following questions

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q10) How long (in kB) is your assembly?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q11) What is the N50?

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q12) How many contigs did abyss generate?

##### Run [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) with your Abyss assembly as the query.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q13) Based on the BLAST result, what organism do you think this genome came from and what kind of genome is it?
- Organism A, genome type A
- Organism B, genome type B
- Organism C, genome type C
- Organism D, genome type D



### BONUS (+2 points)

We did not use all of the options in `velvet` and it can certainly be optimized to produce a better assembly.

Can you change velvet's parameters to yield a genome with a single contig? 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Bonus) What code did you use?

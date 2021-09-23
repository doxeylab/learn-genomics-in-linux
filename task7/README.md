# Task7 - Transcriptomics and detection of differentially expressed genes (DEGs)

In this lab, you will be analyzing RNA-seq data from a study by Aguiar et al. (2020) [here](https://pubmed.ncbi.nlm.nih.gov/31646766/).
This study exposed a human lung cell line (Calu-3 cells) to tobacco smoke, cannabis smoke, and a common drug intervention (LABA/GCS).

You will be measuring and comparing the transcript expression levels between normal untreated cells (controls) and cells exposed to tobacco smoke extract (TSE).
There are 4 TSE samples vs 4 control samples as labeled below.

| Sample ID | Status |
| --------------- | --------------- |
| SRR8451881 | Control |
| SRR8451882 | Control |
| SRR8451883 | Control |
| SRR8451884 | Control |
| SRR8451885 | TSE |
| SRR8451886 | TSE |
| SRR8451887 | TSE |
| SRR8451888 | TSE |

The goal is to identify which genes are up-regulated and down-regulated following tobacco smoke exposure.

### Requirements

#### Command-line tools
* Access to a linux-based OS running BASH
* [Salmon](https://combine-lab.github.io/salmon/)

#### Graphical tools

You will also need to download and install R on your own machine with the following packages

* tximport
* DESeq2 or edgeR


## Getting Started

* Login to your linux environment and create a new folder for your task7

```
mkdir task7  #creates folder
cd task7 #enters into folder
```

## Retrieving the raw data and reference transcriptome

First, download a human reference transcriptome:

```
#download a pre-made reference transcriptome from Gencode
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz
gunzip gencode.v29.transcripts.fa.gz
```

Next, download the RNA-seq data from the public EBI FTP site. We will be downloading 8 samples (each has forward and reverse .fastq reads) - so 16 files total. To download all of these automatically, you can download a list of urls with `wget` as follows:

```
#download the list of urls first - NOTE: THIS STEP CAN TAKE A LONG TIME (~1 hr)
wget https://raw.githubusercontent.com/doxeylab/learn-genomics-in-unix/master/task7/ftp-list.txt

#now create a folder for your data
mkdir data
cd data

#download all of the .fastq files into your data folder
wget -i ../ftp-list.txt

#go back to your task7 folder
cd ..
```


## Transcript quantification with Salmon


Now, before you can measure transcript abundance, you must index your reference transcriptome so that it can be analyzed with `Salmon`

```
salmon index -t gencode.v29.transcripts.fa -i gencode_v29_idx
```


Now, let's measure transcript abundance using `Salmon`. For a single sample with paired-end reads (e.g., `forward_reads.fastq.gz` and `reverse_reads.fastq.gz`, this can be done using the following line:

```
#result will be output to "quants" folder
# -p 6 means that six CPU threads will be used
salmon quant -i gencode_v29_idx -l A -1 forward_reads.fastq.gz -2 reverse_reads.fastq.gz -p 6 -o quants
```

But the above line is just an example for a single sample. Here is a .bash script that will run `Salmon` on all of the 8 samples we have just downloaded.
```
#download bash script
wget https://raw.githubusercontent.com/doxeylab/learn-genomics-in-unix/master/task7/runSalmon.bash

#run bash script. This may take a while...
bash runSalmon.bash

```

## Exploring the transcript counts

Now, the transcript expression levels have been quantified for each of your 8 samples. Look within the `quants` folder and examine the `quant.sf` files that you have produced for each sample.

* Take note of which column contains the transcript id. You can do this by using `head -1` to look at the header of the `quant.sf` file.
* Also take note of which column contains the TPM (transcripts per million) expression level.

Suppose you are interested in the transcript "ENST00000379727.7".

```
#go to your quants/data folder
cd quants

#inspect the expression levels for this transcript
grep "ENST00000379727.7" */quant.sf
```

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q1 - Has this transcript's abundance increased, decreased, or stayed the same following smoke exposure?

Support your answer using statistics. Perform a t-test comparing the expression level of this transcript between the 4 smoke-treated samples versus 4 control samples. Use any program of your choice to do so (R, excel, Google Sheets, etc.).

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q2 -  Is the difference statistically significant (p < 0.01)?

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q3 -  Does this result make sense biologically given existing literature? Provide an answer using 2-3 written sentences, and provide a citation to support your answer.

## Detecting differentially expressed genes (DEGs) in R (Bonus)

Now that you have measured transcript abundance for all samples using `Salmon`, you can perform a differential expression analysis using a tool such as DeSeq2 or edgeR. Here is a rough guide to the steps required:

* Install R on your machine
* Install the `tximport` R package
* Install either the `edgeR` or `deseq2` R package
* Download the quant files produced by Salmon to your local machine
* Following the instructions [here](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Bonus (+2 marks) - Produce a table of the top 10 differentially expressed genes along with their fold-changes and adjusted p-values. Also include the code you used to do so.


---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) above.
Please submit your answers under quizzes on LEARN.




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
mkdir transcriptomics-task  #creates folder
cd transcriptomics-task #enters into folder
```

## Retrieving the raw data and reference transcriptome

First, we need a human reference transcriptome and Salmon index.

This has been done for you already and the files are located at : `/fsys1/data/task4`

If you are curious and would like to know how this was done, see below, but again this is not needed.

```
#download a pre-made reference transcriptome from Gencode
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz
gunzip gencode.v29.transcripts.fa.gz

#index your reference transcriptome so that it can be analyzed with `Salmon`
salmon index -t gencode.v29.transcripts.fa -i gencode_v29_idx


```

Next, we need the the RNA-seq data (8 samples -- fw and rv reads, so 16 total files) from the public EBI FTP site. This has also been downloaded for you. These files were downloaded using the following code:

```
#download the list of urls first - NOTE: THIS STEP CAN TAKE A LONG TIME (~1 hr)
wget https://raw.githubusercontent.com/doxeylab/learn-genomics-in-linux/master/task7/ftp-list.txt

#download all of the .fastq files into your data folder
wget -i ../ftp-list.txt

```


## Transcript quantification with Salmon

Now, you are going to measure transcript abundance using `Salmon`. For a single sample with paired-end reads (e.g., `forward_reads.fastq.gz` and `reverse_reads.fastq.gz`, this could be done using the following line:

```
#result will be output to "quants" folder
# -p 6 means that six CPU threads will be used
salmon quant -i gencode_v29_idx -l A -1 forward_reads.fastq.gz -2 reverse_reads.fastq.gz -p 6 -o quants
```

But the above line is just an example for a single sample. Here is a .bash script that will run `Salmon` on all of the 8 samples we have just downloaded. Run this bash script in your `/transcriptomics-task` folder
```
#download bash script
wget https://raw.githubusercontent.com/doxeylab/learn-genomics-in-linux/master/task7/runSalmon.bash

#run bash script. This may take a few hours...
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

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1 - Has this transcript's abundance increased, decreased, or stayed the same following smoke exposure?

Support your answer using statistics. Perform a t-test comparing the expression level of this transcript between the 4 smoke-treated samples versus 4 control samples. Use any program of your choice to do so (R, excel, Google Sheets, etc.).

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2 -  Is the difference statistically significant (p < 0.05)?


## Detecting differentially expressed genes (DEGs) in R

Now that you have measured transcript abundance for all samples using `Salmon`, you can perform a differential expression analysis using a tool such as DeSeq2 or edgeR. 

On your local machine, install the following R packages:

```
tximport
deseq2
```

Now, on your local machine, open your terminal and download the following quant files produced by Salmon to your local machine

```
scp -r userid@genomics1.private.uwaterloo.ca:~/task4/quants/ .
```

Now, open R and load packages and set working directory:

```
#load required packages
library(tximport)
library(DESeq2)

#go to your folder containing your quant files you just downloaded
setwd("/path/to/quants")

```

Download the gencode reference transcriptome

```
system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.metadata.HGNC.gz")
system("gunzip gencode.v29.metadata.HGNC.gz")
genesymbols = read.delim("gencode.v29.metadata.HGNC")
```

Read the quant files into R

```
files = paste(list.dirs('.', recursive=FALSE),"/","quant.sf",sep='')
#make sure to check your list of files to ensure that this step worked

txi.salmon <- tximport(files, type = "salmon", tx2gene = genesymbols, ignoreAfterBar =T)
```

Run DESeq2 to detect differentially expressed genes between the two categories

```
meta = data.matrix(cbind(files,as.numeric(c(0,0,0,0,1,1,1,1))))
colnames(meta) = c("filenames","category")

dds <- DESeqDataSetFromTximport(txi.salmon, meta, ~as.factor(category))   # no differential design
dds <- DESeq(dds)

res <- results(dds, lfcThreshold=0.5,alpha=0.01)
```

Examine the results

```
summary(res)

```


![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) (2 marks) - Produce a table of the top 10 differentially expressed genes along with their fold-changes and adjusted p-values. Also include the code you used to do so.


---

# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) above.
Please submit your answers under quizzes on LEARN.




# Task2 - Genome Assembly

In this task, you will download raw sequencing data, perform genome assembly, genome annotation, and then compare your assembly to the database using BLAST.

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

Now explore the first few lines of the file using `head`

```
head mt_reads.fastq.gz
```

Now try...

```
less mt_reads.fastq.gz
```

And...

```
tail mt_reads.fastq.gz
```

How many lines are in the file?

```
wc -l mt_reads.fastq.gz
```



## Quality filtering

Before we can assemble a genome, we need to:

* Assess the quality of the sequencing data
* Filter out low-quality reads (this is called quality filtering)

For a quick quality report, you can use the program `fastqc`.
The command below will analyze the mt_reads.fastq file and place the results in the 'qcresults' folder.

```
fastqc -o qcresults mt_reads.fastq
```

Now open the HTML file that was generated in firefox. The ampersand (&) will run the command in the background so you can still use the terminal.

```
firefox qcresults/mt_reads_fastqc/fastqc_report.html &
```

###Splitting the barcodes



## Genome Assembly

...

## Genome Annotation

...

## BLAST




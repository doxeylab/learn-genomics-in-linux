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

gunzip mt_read.fastq.gz

```

Now explore the first few lines of the file using `head`

```
head mt_reads.fastq.gz
```

Now try

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

Before we can assemble a genome, we need to filter out low quality reads.


## Genome Assembly

...

## Genome Annotation

...

## BLAST




# Task8 - Computational metagenomics tutorial (16S amplicon and shotgun metagenomics)

In this lab, you will be analyzing 16S and metagenomic data from a study by Lobb et al. (2020) [here](https://pubmed.ncbi.nlm.nih.gov/32345738/).
This study examined the microbial communities of decomposing fish in local rivers near Waterloo, ON, Canada.

There are 52 samples with the following metadata. 

| Sample ID | Name |
| --------------- | --------------- |
| SRS6112281 | WW1.1 |
| SRS6112282 | WW1.2 |
| SRS6112283 | WW1.3 |
| SRS6112267 | EE1.1 |
| SRS6112268 | EE1.2 |
| SRS6112279 | EE1.3 |
| SRS6112285 | WW4.1 |
| SRS6112284 | WW4.2 |
| SRS6112286 | WW4.3 |
| SRS6112293 | WE4.1 |
| SRS6112295 | WE4.2 |
| SRS6112296 | WE4.3 |
| SRS6112290 | EE4.1 |
| SRS6112302 | EE4.2 |
| SRS6112304 | EE4.3 |
| SRS6112271 | EW4.1 |
| SRS6112272 | EW4.2 |
| SRS6112273 | EW4.3 |
| SRS6112287 | WW8.1 |
| SRS6112288 | WW8.2 |
| SRS6112289 | WW8.3 |
| SRS6112297 | WE8.1 |
| SRS6112298 | WE8.2 |
| SRS6112299 | WE8.3 |
| SRS6112305 | EE8.1 |
| SRS6112306 | EE8.2 |
| SRS6112307 | EE8.3 |
| SRS6112274 | EW8.1 |
| SRS6112275 | EW8.2 |
| SRS6112276 | EW8.3 |
| SRS6112291 | WW10.1 |
| SRS6112292 | WW10.2 |
| SRS6112294 | WW10.3 |
| SRS6112300 | WE10.1 |
| SRS6112301 | WE10.2 |
| SRS6112303 | WE10.3 |
| SRS6112308 | EE10.1 |
| SRS6112269 | EE10.2 |
| SRS6112270 | EE10.3 |
| SRS6112277 | EW10.1 |
| SRS6112278 | EW10.2 |
| SRS6112280 | EW10.3 |

Our goal will be to perform taxonomic profiling by analyzing the 16S dataset. We will then perform targeted metagenome assembly to investigate a sample of interest.

### Requirements

#### Command-line tools
* Access to a linux-based OS running BASH
* kraken2 and Bracken
* metaspades

#### Graphical tools

You will also need to download and install R on your own machine with the following packages

* ggplot2
* pheatmap
* reshape2
* viridisLite


## Getting Started

* Login to your linux environment and create a new folder for your task7

```
mkdir task8  #creates folder
cd task8 #enters into folder
```

## Retrieving the raw data

The data has already been downloaded for you, and is located in the `/fsys/data/lobb-et-al/` folder

If you're curious, the original data was downloaded from the NCBI SRA using this command:
```
fastq-dump --split-files SRS6112303 SRS6112301 SRS6112300 SRS6112299 SRS6112298 SRS6112297 SRS6112296 SRS6112295 SRS6112293 SRS6112294 SRS6112292 SRS6112291 SRS6112289 SRS6112288 SRS6112287 SRS6112286 SRS6112284 SRS6112285 SRS6112283 SRS6112282 SRS6112281 SRS6112280 SRS6112278 SRS6112277 SRS6112276 SRS6112275 SRS6112274 SRS6112273 SRS6112272 SRS6112271 SRS6112270 SRS6112269 SRS6112308 SRS6112307 SRS6112306 SRS6112305 SRS6112304 SRS6112302 SRS6112290 SRS6112279 SRS6112268 SRS6112267 SRS6098991 SRS6098990 SRS6098989 SRS6098988 SRS6098999 SRS6098998 SRS6098997 SRS6098996 SRS6098995 SRS6098994 SRS6098993 SRS6098992 SRS6098987 SRS6098986
```


## Taxonomic classification of 16S reads using Kraken2

Tools such as `QIIME2` and `Mothur` are common for analyzing 16S rRNA sequences. For this tutorial, we will be using a different tool called `Kraken2`.

Suppose we wanted to analyze a single sample (e.g., SRS6112303). We can do so with the following Kraken2 command:

```
CLASSIFICATION_LVL=G  # this will set an environmental variable for the taxonomic level of classification desired (G = "Genus", S = "Species", etc.)
krakenDB=/data/krakendb/16S_Greengenes_k2db/  #this is the location of the kraken2 database you want to use for classification
fastq1=/fsys1/data/lobb-et-al/SRS6112303_1.fastq
fastq2=/fsys1/data/lobb-et-al/SRS6112303_2.fastq

kraken2 --db $krakenDB --paired --report report.txt --output kraken.out $fastq1 $fastq2

bracken -d $krakenDB -l $CLASSIFICATION_LVL -i report.txt -o bracken.out
```

`bracken.out` will look like this:

```
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Parabacteroides distasonis	2601	S	8	53	61	0.00554
Parabacteroides gordonii	2602	S	1	37	38	0.00347
Bacteroides fragilis	2596	S	11	2044	2055	0.18386
Bacteroides uniformis	2599	S	1	9	10	0.00094
Clostridium pasteurianum	2769	S	278	346	624	0.05586
Clostridium perfringens	2770	S	146	986	1132	0.10128
Clostridium subterminale	2772	S	78	2011	2089	0.18691
Clostridium bowmanii	2764	S	13	64	77	0.00692
Clostridium butyricum	2765	S	10	102	112	0.01009
Clostridium neonatale	2768	S	3	33	36	0.00328
Alkaliphilus transvaalensis	2762	S	3	1	4	0.00037
Sporomusa polytropa	2800	S	4	98	102	0.00913
Veillonella dispar	2801	S	22	255	277	0.02478
...
...

```

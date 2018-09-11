# Task1b - Command-line BLAST

In this task, you will learn how to run BLAST in the command-line. You will download a genome and proteome and search them for a gene and protein of interest, respectively.

### Requirements

* Access to a linux-based OS running BASH
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

---

## Getting Started

Login to your linux environment as you did in task1. e.g.,

```
ssh -i /path/to/your/.ssh/publickey yourUserName@remoteIP
```

Create a new folder for this task

```
mkdir blastTask  #creates folder
cd blastTask #enters into folder
```

## Retrieving the raw data

There are several ways to download data. Two common tools are `curl` and `wget`.

The following command will download a genome (DNA sequence data) and proteome (translated AA sequence data) from the NCBI.

```
wget ...
wget ...

```

## Formatting the genome and proteome for BLAST

BLAST requires that the genome and proteome be formatted. What BLAST is doing here is 'indexing' the sequence data, which chops sequences into their constitutent k-mer fragments and stores this data to facilitate rapid database searching.

To format the data, it must be in FASTA format, which looks like this:

```
>sequenceName
ATCATCGTACGATCGATCATCG
ATCATCTAGCTATGCATCGTAC
AATGCACGTACGTATCGACTCT
```
or for protein
```
>sequenceName
MARSTLEDQIIAAIDL
etc. etc.
```

Take a look at your input files to ensure that they adhere to this format using `head` or `less`

To format the genome for BLAST, do the following

```
makeblastdb -in <genome.fna> -dbtype 'nucl'
```

To format the proteome for BLAST, do the same thing but change the 'dbtype'

```
makeblastdb -in <proteome.fna. -dbtype 'prot'
```


# ASSIGNMENT QUESTIONS

Your assignment will be to write a series of shell commands to answer the following questions. Please submit the code you used as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.

* Pick a protein of interest (POI) (from any organism)

* Now pick a genome from a different organism (these can be found in the the NCBI ftp directory) where you suspect to find a homolog of your protein

* Download the genome (.fna) and predicted proteome (.faa) for the organism above

* Search the genome for your POI using `tblastn` and the proteome for your POI using `blastp`

* Paste your top result for each.

* Have you a identified true homolog of your POI? Justify your answer using the E-value and sequence identity.

* How many homologs did you detect? Again, justify your answer by pasting the relevant output from your BLAST result.




#### Congratulations. You are now finished Task 1b.


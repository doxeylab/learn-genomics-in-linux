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

Create a new project folder for this task

```
mkdir blastTask  #creates folder
cd blastTask #enters into folder
```

## Retrieving the raw data

There are several ways to download data. Two common tools are `curl` and `wget`.

The following command will download a genome (DNA sequence data) and proteome (translated AA sequence data) from the NCBI.
The NCBI houses its genomic data within an FTP directory - [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank)
It is further subdivided - e.g., bacterial genomes can be found [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria)

We will be working with the genome of Prochlorococcus marinus, which is an abundant marine microbe and possibly the most abundant genus on earth. First, explore its FTP directory [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1)

Within this folder, there are a number of files including a:

* ...genomic.fna.gz file -> this will uncompress into a .fna (fasta nucleic acid) file
* ...protein.faa.gz -> this will uncompress into a .faa (fasta amino acid) file

It should be clear what these two files contain.

There is also another file called:

* ...genomic.gbk.gz -> this will uncompress into a .gbk file (genbank file)


Download these files, uncompress them, and explore them (with `less` for example).

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_genomic.fna.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_genomic.gff.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_protein.faa.gz

```

<b>Q1) Look at the .gbk file. What information does this contain? </b>


## Formatting the genome and proteome for BLAST

Next, we are going to do a BLAST search against the genome (.fna) and proteome (.faa) that you have downloaded.

In this case, the genome and proteome will be the <b>DATABASE</b> you are searching against.
The <b>QUERY</b> sequences (a gene and a protein) can be anything you like, but here are some suggestions:

* e.g. query gene sequence: E. coli 16S ribosomal RNA - copy and paste this into a new text file
``
>J01859.1 Escherichia coli 16S ribosomal RNA, complete sequence
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT
AACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATG
GAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCG
GGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACG
ATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGG
CAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT
CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCG
CAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAAT
TACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAAC
TGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGT
AGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCG
TGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCC
TTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACT
CAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCT
TACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGC
TGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCT
TTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGA
CGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGA
CCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATG
AAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCG
CCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTT
TGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
A
``

* e.g. query protein sequence: 

``curl https://www.uniprot.org/uniprot/B7LA79.fasta > e.coli.l7.faa``


When doing a BLAST, the query can be in FASTA format, but the database needs to be <b>formatted</b> for BLAST. This is done with the `makeblastdb` command. This tool sets up an 'indexed' database for BLAST, which chops sequences into their constitutent k-mer fragments and stores this data to facilitate rapid database searching.

Let's set up a BLAST database for the proteome

```
makeblastdb -in GCA_000007925.1_ASM792v1_protein.faa -dbtype 'prot'
```

... and now the genome as well

```
makeblastdb -in GCA_000007925.1_ASM792v1_genomic.faa -dbtype 'nucl'
```


You can see that the `-dbtype` parameter defines whether the input FASTA file is for protein or nucleotide sequences.



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


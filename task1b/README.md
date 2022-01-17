# Task1b - Command-line BLAST

In this task, you will learn how to run BLAST in the command-line. You will download a genome and proteome and search them for a gene and protein of interest, respectively.

### Requirements

* Access to a linux-based OS running BASH
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

---

## Getting Started

Login to your linux environment as you did in task1, either through the browser or via `ssh`.

Create a new project folder for this task

```
mkdir blastTask  #creates folder
cd blastTask #enters into folder
```

## Retrieving genomic data from the NCBI

There are several ways to download data. Two common tools are `curl` and `wget`.
You can also simply copy and paste sequence data into a file using `nano` or `pico` or other command-line text-editors. More advanced ones are `vim` and `emacs`.

The following exercise will download a genome (DNA sequence data) and proteome (translated AA sequence data) from the NCBI.
The NCBI houses its genomic data within an FTP directory - [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/)

We will be working with the genome of [<i>Prochlorococcus marinus</i>](https://en.wikipedia.org/wiki/Prochlorococcus), which is an abundant marine microbe and possibly the most abundant bacterial genus on earth. First, explore its FTP directory [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/925/GCA_000007925.1_ASM792v1/)

Within this folder, there are a number of files. It is important that you familiarize yourself with these files and their contents.
Files within a genbank genomic ftp directory include:

* ...genomic.fna.gz file -> this will uncompress into a .fna (fasta nucleic acid) file
* ...protein.faa.gz -> .faa (fasta amino acid) file

It should be clear what these two files contain.

There is also another file called:

* ...genomic.gff.gz -> .gff file (generic feature format)


Download these files, uncompress them, and explore them (with `less` for example).

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/925/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_genomic.fna.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/925/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_genomic.gff.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/925/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_protein.faa.gz

```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1 - Examine the contents of the three files you have just downloaded. What information does each contain?


## Command-line BLAST

### Setting up your query sequence(s)

Next, you are going to do a BLAST search against the genome (.fna) and proteome (.faa) that you have downloaded.

In this case, the genome and proteome will set up as BLAST <b>DATABASES</b> that you are searching against.
The BLAST <b>QUERY</b> sequences (a gene and a protein) can be anything you like (examples below).

* e.g. a query gene sequence: E. coli 16S ribosomal RNA - copy and paste this into a new text file

```
>J01859.1 Escherichia coli 16S ribosomal RNA, complete sequence
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA
```

* e.g. query protein sequence: 

```
>sp|B7LA79|RL7_ECO55 50S ribosomal protein L7/L12 OS=Escherichia coli (strain 55989 / EAEC) OX=585055 GN=rplL PE=3 SV=1
MSITKDQIIEAVAAMSVMDVVELISAMEEKFGVSAAAAVAVAAGPVEAAEEKTEFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK
```

Note: here is a quick way to download the above query protein sequence from Uniprot and rename it in one command

```
curl https://www.uniprot.org/uniprot/B7LA79.fasta > e.coli.l7.faa
```


### Formatting the genome and proteome for BLAST

When doing a BLAST search, the query can be in FASTA format, but the database needs to be <b>formatted</b> for BLAST. This is done with BLAST's `makeblastdb` command. This tool sets up an 'indexed' database for BLAST, which chops sequences into their constitutent k-mer fragments and stores this data to facilitate rapid database searching.

Let's set up a BLAST database for the proteome

```
makeblastdb -in GCA_000007925.1_ASM792v1_protein.faa -dbtype 'prot'
```

... and now the genome as well

```
makeblastdb -in GCA_000007925.1_ASM792v1_genomic.fna -dbtype 'nucl'
```

You can see that the `-dbtype` parameter defines whether the input FASTA file is for protein or nucleotide sequences.

`makeblastdb` and other BLAST tools have lots of additional parameters as well that can be customized based on a users needs.
Let's explore some more.

```
# to look at command usage and parameter options
makeblastdb -help
```

#### Advanced: retrieving specific entries and regions from your BLAST database

One of the useful parameters here is the `-parse_seq_ids` flag. If this option is set, this makes it very easy to retrieve specific sequences from the database using their name or id. e.g.,

```
makeblastdb -in GCA_000007925.1_ASM792v1_genomic.fna -dbtype 'nucl' -parse_seqids
makeblastdb -in GCA_000007925.1_ASM792v1_protein.faa -dbtype 'prot' -parse_seqids
```

And now if you want to print out to the screen the sequence of protein `AAP99047.1`, you can use the `blastdbcmd` program like this:

```
blastdbcmd -entry AAP99047.1 -db GCA_000007925.1_ASM792v1_protein.faa
```
This will output:
>\>AAP99047.1 DNA polymerase III beta subunit [Prochlorococcus marinus subsp. marinus str. CCMP1375]
MKLVCSQIELNTALQLVSRAVATRPSHPVLANVLLTADAGTGKLSLTGFDLNLGIQTSLSASIESSGAITVPSKLFGEII
SKLSSESSITLSTDDSSEQVNLKSKSGNYQVRAMSADDFPDLPMVENGAFLKVNANSFAVSLKSTLFASSTDEAKQILTG
VNLCFEGNSLKSAATDGHRLAVLDLQNVIASETNPEINNLSEKLEVTLPSRSLRELERFLSGCKSDSEISCFYDQGQFVF
ISSGQIITTRTLDGNYPNYNQLIPDQFSNQLVLDKKYFIAALERIAVLAEQHNNVVKISTNKELQILNISADAQDLGSGS
ESIPIKYDSEDIQIAFNSRYLLEGLKIIETNTILLKFNAPTTPAIFTPNDETNFVYLVMPVQIRS

If you want a specific region (e.g., the first 10 amino acids) from this entry, you can use the `-range` parameter

```
blastdbcmd -entry AAP99047.1 -db GCA_000007925.1_ASM792v1_protein.faa -range 1-10
```
This will output:
>\>AAP99047.1 DNA polymerase III beta subunit [Prochlorococcus marinus subsp. marinus str. CCMP1375]
MKLVCSQIEL


### Performing a BLAST search

There are several different flavors of BLAST. Each is run as a separate command:

* `blastp` - protein query vs protein database
* `blastn` - nucleotide query vs nucleotide database
* `blastx` - nucleotide query (translated) vs protein database
* `tblastn` - protein query vs nucleotide (translated) database
* `tblastx` - nucleotide query (translated) vs nucleotide database (translated)

To run a `blastp` search using the protein query (defined by `-query` parameter) and protein database (defined by `-db` parameter) you have set up, do the following:

```
blastp -query e.coli.l7.faa -db GCA_000007925.1_ASM792v1_protein.faa
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2 - How many significant (E < 0.001) hits did you get?

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q3 - What is the sequence identity (percentage) of your top BLAST hit?

* Repeat the same BLAST search you did for Q2 but using the genomic sequence as the search database.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q4 - Compare your result to the previous search. Which of the following statements is most correct:

```
* There were no significant BLAST matches.
* BLAST detected the same protein as the top hit. However, the alignment was shorter.
* BLAST detected the same protein as the top hit. However, the alignment was not significant.
* BLAST detected a different protein as the top hit.
```


* Suppose you have sequenced the following fragment of DNA:

```
ACTGGCATTGATAGAACAACCATTTATTCGAGATAGTTCAATTACTGTAGAGCAAGTTGTAAAACA
```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q5 - Search for this fragment of DNA in the genome of Prochlorococcus marinus subsp. marinus str. CCMP1375. What did you find?

```
* I found an exact match to the sequence.
* I did not find a good match to the sequence.
* I found a good match to the sequence with 1 mutation.
* I found a good match to the sequence with 2 mutations.

```

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q6 - Based on your BLAST result, what is the likely function of this fragment of DNA?

```
* It is impossible to say.
* It is part of a gene encoding Translation elongation factor Ts.
* It is a segment of a protein.
* It is a non-coding sequence.
```


# ASSIGNMENT QUESTIONS

<i>* Complete questions 1-6 above and submit your answers on LEARN.</i>


#### Congratulations. You are now finished Task 1b.


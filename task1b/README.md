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
You can also simply copy and paste sequence data into a file using `nano` or `pico` or other unix text-editors. More advanced ones are `vim` and `emacs`.

The following command will download a genome (DNA sequence data) and proteome (translated AA sequence data) from the NCBI.
The NCBI houses its genomic data within an FTP directory - [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank)
It is further subdivided - e.g., bacterial genomes can be found [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria)

We will be working with the genome of [<i>Prochlorococcus marinus</i>](https://en.wikipedia.org/wiki/Prochlorococcus), which is an abundant marine microbe and possibly the most abundant bacterial genus on earth. First, explore its FTP directory [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1)

Within this folder, there are a number of files. It is important that you familiarize yourself with these files and their contents.
Files within a genbank genomic ftp directory include:

* ...genomic.fna.gz file -> this will uncompress into a .fna (fasta nucleic acid) file
* ...protein.faa.gz -> .faa (fasta amino acid) file

It should be clear what these two files contain.

There is also another file called:

* ...genomic.gbk.gz -> .gbk file (genbank file)


Download these files, uncompress them, and explore them (with `less` for example).

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_genomic.fna.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_genomic.gff.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Prochlorococcus_marinus/reference/GCA_000007925.1_ASM792v1/GCA_000007925.1_ASM792v1_protein.faa.gz

```

:large_blue_circle: <b>1 - Look at the .gbk file. What information does this contain? </b>


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
makeblastdb -in GCA_000007925.1_ASM792v1_genomic.faa -dbtype 'nucl'
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

:large_blue_circle: <b>2 - How many significant (E < 0.001) hits did you get?</b>

:large_blue_circle: <b>3 - Copy and paste the alignment for the top hit.</b>

:large_blue_circle: <b>4 - Search this protein against the genomic sequence. Do you get the same result? Describe your answer.</b>


# ASSIGNMENT QUESTIONS

:large_blue_circle: * Answer questions 1 - 4 above

:large_blue_circle: 5 - Pick a protein of interest (POI) (from any organism)

:large_blue_circle: 6 - Now pick a genome from a different organism (these can be found in the the NCBI ftp directory) where you suspect to find a homolog of your protein

:large_blue_circle: 7 - Download the genome (.fna) and predicted proteome (.faa) for the organism above

:large_blue_circle: 8 - Search the genome for your POI using `tblastn` and the proteome for your POI using `blastp`

:large_blue_circle: 9 - Paste your top result for each.

:large_blue_circle: 10 - Have you a identified true homolog of your POI? Justify your answer using the E-value and sequence identity.

:large_blue_circle: 11 - How many homologs did you detect? Again, justify your answer by pasting the relevant output from your BLAST result.


Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.


#### Congratulations. You are now finished Task 1b.


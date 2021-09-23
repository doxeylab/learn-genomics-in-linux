# Task4 - Synteny comparison of genomes

This task is a tutorial on structural comparison of genomes using synteny mapping.

### Requirements

* Access to a linux-based OS running BASH
* [BLAST](http://blast.ncbi.nlm.nih.gov/)
* [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) *
* [Mauve](http://darlinglab.org/mauve/download.html) (optional) *

## Installation

If you do not already have access to a GUI running the graphical software listed above (*), please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

All software used are available for Mac/Windows/Linux.

---

## Getting Started

* Login to your linux environment and create a new folder for task4

```
mkdir task4  #creates folder
cd task4 #enters into folder
```

## Retrieving the raw data

Copy the genome from task 2 you assembled with `abyss`

```
cp ../task2/abyss-assembly-contigs.fa . 
```

You will be comparing this genome to another related genome from L. terrestris. Download this genome.

```
wget https://github.com/doxeylab/learn-genomics-in-unix/raw/master/task4/l-terrestris.genome.fa
```

Make BLAST databases for both.

```
makeblastdb -in abyss-assembly-contigs.fa -dbtype nucl
makeblastdb -in l-terrestris.genome.fa -dbtype nucl
```

Now BLAST one genome against the other with the following command. Note that you are using BLAST's `-outfmt 6` parameter which outputs the BLAST result as a table (which you are writing to `blastresults.tab`). You will be using this table to visualize the synteny between these two genomes.

```
blastn -outfmt 6 -db abyss-assembly-contigs.fa -query l-terrestris.genome.fa >blastresults.tab
```

Now, download to your local machine the following files:

* abyss-assembly-contigs.fa
* l-terrestris.genome.fa
* blastresults.tab

Open the `act` program that is packaged with `artemis` and input these three files.

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q1) Paste a screenshot of your result.

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q2) Describe the synteny patterns that you are observing. Is there 1:1 synteny involving the entire genome or have there been genomic rearrangements?

---

## Working with your own dataset

Next, find two related genomes (i.e., different strains of same species)  from the [NCBI FTP directory](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/)

* Repeat the analyses above to perform a structural genome comparison.

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q3) Paste a screenshot of your result.

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Q4) Describe the synteny patterns that you are observing.


## Multiple genome alignment with MAUVE -- Bonus (+1)

This is for bonus marks.

Want to try aligning/comparing more than two genomes? 

* Download/install [mauve](http://darlinglab.org/mauve/download.html) to your local machine

* Select three or more genomes of interest

* Open the sequences in `Mauve` and align them

* Visualize the multiple alignment

![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) Bonus) Paste a screenshot of your result.


---


# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![question](https://github.com/doxeylab/learn-genomics-in-unix/blob/master/questionbox.png) above.
Please submit the code you used (when required) as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.


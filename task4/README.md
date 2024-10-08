# Task4 - Synteny comparison of genomes

This task is a tutorial on structural comparison of genomes using synteny mapping.

### Requirements

* Access to a linux-based OS running BASH
* [BLAST](http://blast.ncbi.nlm.nih.gov/)
* [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) * download this graphical software onto your own machine. Again, see [here](installingArtemis.md) for further instructions
* [Mauve](http://darlinglab.org/mauve/download.html) (optional) * download this graphical software onto your own machine

## Installation

<!--If you do not already have access to a GUI running the graphical software listed above (*), please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.-->Please install the graphical software on your local machine.

All software used are available for Mac/Windows/Linux.

---

## Getting Started

* Login to your linux environment and create a new folder for task4.

```
mkdir task4  #creates folder
cd task4 #enters into folder
```

## Retrieving the raw data

* Copy the genome from task2 you assembled with `abyss`

```
cp ../task2/abyss-assembly-contigs.fa . 
```

* You will be comparing this genome to another related genome from <i>L. terrestris</i>. Download this genome.

```
wget https://github.com/doxeylab/learn-genomics-in-linux/raw/master/task4/l-terrestris.genome.fa
```

* Make BLAST databases for both.

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

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q1) Paste a screenshot of your result.  (3 marks)

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q2) Describe the synteny pattern that you are observing. Do you think genomic rearrangements have taken place or is there a strong pattern of shared synteny between both genomes? (2 marks) See [shared synteny](https://en.wikipedia.org/wiki/Synteny#Shared_synteny).

To help you with this question, consider two genome sequences composed of four genes A-D. One genome has gene order A,B,C,D and the second genome has gene order A,C,B,D. There has clearly been a genomic rearrangement here because C and B have switched places.

But now suppose the genomes are (A,B,C,D) and (C,D,A,B). If these are linear chromosomes, then a rearrangement has taken place, but what if they are circular?

And lastly, now suppose we compare (A,B,C,D) to its reverse complement which will appear to be in the order (D,C,B,A). This may look like an inversion in artemis, but one of the two strands just needs to be flipped so that we are comparing the genomes in the same orientation.

---

## Working with your own dataset

Next, find two related genomes (e.g., different strains of same species)  from the [NCBI Genome Database](https://www.ncbi.nlm.nih.gov/datasets/genome/). <!--[NCBI FTP directory](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/)-->

* Repeat the analyses above to perform a structural genome comparison.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q3) Paste a screenshot of your result. (3 marks)

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Q4) Describe the synteny patterns that you are observing. (2 marks)


## Multiple genome alignment with Mauve -- Bonus (+1)

This is for bonus marks.

Want to try aligning/comparing more than two genomes? 

* Download/install [Mauve](http://darlinglab.org/mauve/download.html) to your local machine.

* Select three or more genomes of interest.

* Open the sequences in `Mauve` and align them.

* Visualize the multiple alignment.

![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) Bonus) Paste a screenshot of your result.


---


# ASSIGNMENT QUESTIONS

The questions for this task are indicated by the lines starting with ![question](https://github.com/doxeylab/learn-genomics-in-linux/raw/master/questionbox.png) above.
Please submit the code you used (when required) as well as the answers to the questions. Submit your assignment to a dropbox on LEARN as a .docx, .txt, or .pdf file.


# Task4 - Synteny comparison of genomes

This task is a tutorial on structural comparison of genomes using synteny mapping.

### Requirements

* Access to a linux-based OS running BASH
* [BLAST](http://blast.ncbi.nlm.nih.gov/)
* [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/)
* [Mauve](http://darlinglab.org/mauve/download.html)

## Installation

If you do not already have access to a GUI running the graphical software listed above, please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

All software used are available for Mac/Windows/Linux.

---

## Getting Started

* Login to your linux environment as you did in task1.

* Create a new folder for your task4

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





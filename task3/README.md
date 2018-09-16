# Task3 - Genome Annotation

In this task, you will perform more detailed annotation of a genome and explore different types of ontologies.

### Requirements

* Access to a linux-based OS running BASH
* Prokka

## Installation

If you do not already have access to a GUI running the graphical software listed above, please install the software on your local machine. Once locally installed, you can download results off the linux server and locally visualize them on your own system.

All software used are available for Mac/Windows/Linux.

---

## Getting Started

Login to your linux environment as you did in task1. e.g.,

```
ssh -i /path/to/your/.ssh/publickey yourUserName@remoteIP
```

Create a new folder for your task3

```
mkdir task3  #creates folder
cd task3 #enters into folder
```

## Retrieving the raw data

Copy the genome you assembled from task2

```
cp ../task2/abyss-assembly-contigs.fa . 
```

## Genome Annotation I

By marking the ORFs in your genome (given a min size threshold), you have essentially performed a simple gene finding algorithm. However, there are more advanced ways of gene-finding that take additional criteria into account.

A popular genome annotation tool for prokaryotic genomes is [`prokka`](https://github.com/tseemann/prokka).
`prokka` automates a series of genome annotation tools and is simple to run. It has been installed for you on the server.

Type

```
prokka abyss-assembly-contigs.fa
```

* Now, download the .gbk file that was produced and view it in Artemis

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q1) Are the predicted gene locations consistent with your earlier ORF predictions from task 2? 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q2) Why are there vertical black lines in the middle of predicted ORFs?

* Quit artemis and change your artemis 'Options' to better reflect the source of this genome. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) Q3) Does this correct the issues above? Paste a screenshot.

## Genome Annotation II

Next, let's perform genome annotation on a larger scale.

Download the E. coli K12 genome from task1

```
wget
```



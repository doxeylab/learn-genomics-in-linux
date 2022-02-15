# Task8 - Analysis of 16S amplicon sequencing data using Kraken2/Bracken

In this lab, you will be analyzing 16S and metagenomic data from a study by Lobb et al. (2020) [here](https://pubmed.ncbi.nlm.nih.gov/32345738/).
This study examined the microbial communities of decomposing fish in local rivers near Waterloo, ON, Canada.

There are 48 16S rRNA samples with the following metadata. 

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

Our goal will be to perform taxonomic profiling of these 16S rRNA datasets using Kraken2/Bracken.

### Requirements

#### Command-line tools
* Access to a linux-based OS running BASH
* kraken2 and Bracken
* metaspades

#### Graphical tools

You will also need to download and install [R](https://www.r-project.org) on your own machine with the following packages

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


## Quality filtering

We have previously covered the use of tools such as `fastqc` and `trimmomatic` to quality filter our dataset. QC is a required step for any high-throughput sequencing pipeline, but for simplicity we will skip it for the purposes of this tutorial.

## Taxonomic classification of 16S reads using Kraken2

### Analyzing single samples

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

### Analyzing many samples

But these commands run Kraken2/Bracken on only a single sample. What do we do if we want to run them on all samples?

First, you need to have a list of the samples you want to analyze. This has been done for you with the file at `/fsys1/data/lobb-et-al/files.txt`.

Then, we will create a bash script called `runAll.bash` with the following contents.

```
#!/bin/bash

# $1 is the file containing the list of samples
# $2 is the classification level

CLASSIFICATION_LVL=$2

while IFS=$'\t' read sample
do 
    echo "processing sample $sample"

    kraken2 --db /data/krakendb/16S_Greengenes_k2db/ --paired --report $sample.$CLASSIFICATION_LVL.kraken --output $sample.$CLASSIFICATION_LVL.kraken.out /fsys1/data/lobb-et-al/${sample}_1.fastq /fsys1/data/lobb-et-al/${sample}_2.fastq

    bracken -d /data/krakendb/16S_Greengenes_k2db -l $CLASSIFICATION_LVL -i $sample.$CLASSIFICATION_LVL.kraken -o $sample.$CLASSIFICATION_LVL.bracken.out

done < $1
```

Let's now run `runAll.bash` to apply Kraken2 and Bracken to all of our 16S samples.

```
#first let's create a new folder
mkdir order_classification

cd order_classification

# this will perform taxonomic classification at the Order level
bash runAll.bash /fsys1/data/lobb-et-al/files.txt O  # now wait a while....

```

Once completed, you will see that your folder is full of .bracken output files.
To merge these together into a single file containing bracken output for all your samples, do the following:


```
python2.7 /usr/local/bin/Bracken-2.5/analysis_scripts/combine_bracken_outputs.py --files $(ls *.bracken.out) -o combined.order.out
```


## Plotting results in R

Now, download your `combined.order.out` file to your local computer and load [R](https://www.r-project.org/) for further analysis and plotting.

First load these libraries. Install these first if they are not already installed.

```
library(ggplot2)
library(reshape2)
library(viridisLite)
```

Next, load your data

```
tb = read.delim("combined.order.out",header=T,row.names=1)

#Bracken output has _frac and _num columns. We will just analyze the _num columns.
tbp = tb[,grep("_num",colnames(tb))]

#Transpose the table
tbp <- t(tbp)

#Convert to proportions
tb_prop<-as.data.frame(round(prop.table(as.matrix(tbp), 1) * 100,1))

#Choose a selection of taxa with a % > 3 (Note: might have to play around with this until you get a reasonable number of taxa to display)
tb_sub <- tb_prop[,apply(tb_prop, 2, function(x) max(x, na.rm = TRUE))>3]

```

For plotting, we have to do a few more modifications to the data matrix
```
#Melt the dataframe for plotting
tbm <- as.data.frame(melt(as.matrix(tb_sub)))

#fix labels
tbm[,1] = within(tbm, Var1<-data.frame(do.call('rbind', strsplit(as.character(Var1), '.', fixed=TRUE))))[,1][,1]

#Turn 0s into NAs
tbm[tbm == 0] <- NA

#Set the order of the taxa on the plot (Note: optional)
tbm$Var2 <- factor(tbm$Var2, levels = row.names(as.table(sort(colMeans(tb_sub)))))

```

Now, we can plot using `ggplot2`. Note: the following ggplot command is very parameter-rich, and it can be a lot simpler than this.

```
ggplot(tbm, aes(Var1,Var2,size = value,fill=value), colsep=c(1:100), rowsep=(1:100), sepwidth=c(5,1)) + geom_point(shape = 21, alpha=0.4) + ggtitle("") + xlab("") + ylab("") + theme(axis.text = element_text(colour= "black", size = 12), text = element_text(size=15), axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))+ scale_size_area(max_size = 15,guide="none") + labs(fill="Relative\nfrequency (%)") + scale_fill_viridis_c()
```
This should produce the following plot:

![](https://github.com/doxeylab/learn-genomics-in-linux/blob/master/task8/bubbleplot.png)

We can also create a barplot by doing the following:

```
ggplot(tbm, aes(fill=Var2, y=value, x=Var1)) + 
    geom_bar(position="fill", stat="identity", col="grey50") +
	scale_y_continuous(labels=scales::percent) +
	xlab("") + ylab("Relative frequency") + labs(fill="Family") +
	theme(axis.text.x= element_text(angle = 90, hjust = 1))

```

... which should produce:

![](https://github.com/doxeylab/learn-genomics-in-linux/blob/master/task8/barplot1.png)


## Adding in metadata annotations

The last plot contains sample names, but let's replace these names with annotations from the metadata that are more informative.

First, make sure you have a `metadata.txt` text file that contains the 48 samples (column 1) and names (column 2) listed at the top of the page.
It should look like this:

```
SRS6112281 WW1.1
SRS6112282 WW1.2
SRS6112283 WW1.3
SRS6112267 EE1.1
SRS6112268 EE1.2
SRS6112279 EE1.3
...
```

Now, load in your metadata
```
metadata = read.table("metadata.txt")
```

Now, let's subset our data matrix to include only the metadata samples, and let's also re-order the variables so that they plot in the desired order

```
#subset tbm to only include those samples in metadata
tbm = tbm[which(tbm[,1] %in% metadata[,1]),]
tbm[,1] = metadata[match(tbm[,1],metadata[,1]),2]
tbm$Var1 = factor(tbm$Var1,levels= metadata[,2])

ggplot(tbm, aes(fill=Var2, y=value, x=Var1)) + 
    geom_bar(position="fill", stat="identity", col="grey50") +
	scale_y_continuous(labels=scales::percent) +
	xlab("") + ylab("Relative frequency") + labs(fill="Family") +
	theme(axis.text.x= element_text(angle = 90, hjust = 1))
```

This should produce:

![](https://github.com/doxeylab/learn-genomics-in-linux/blob/master/task8/barplot2.png)


How does this result compare to the result from Lobb et al. (2020) [here](https://pubmed.ncbi.nlm.nih.gov/32345738/) ?


Lastly, let's create a heatmap and add in an annotation category

```
library(pheatmap)

# convert the tbm table back to a 2D matrix
tb = acast(tbm, Var1 ~ Var2,value.var='value',fill=0)
tb = t(tb)  #transpose

# let's split the names (EE, WW, etc.) into a matrix that we can use as annotations
annot = data.frame(do.call("rbind", strsplit(as.character(metadata[,2]), "", fixed = TRUE)))[,c(1,2)]
rownames(annot) = metadata[,2]
colnames(annot) = c("Fish_Origin","Water_Origin")

# specify the colors
ann_colors = list(
    Fish_Origin = c(W = "#EBEBEB", E = "#424242"),
    Water_Origin = c(W = "#EBEBEB", E = "#424242")
)

# plot
pheatmap(tb,annotation_col=annot,cluster_cols=F,annotation_colors=ann_colors,color = viridis(1000))
```

This should produce:

![](https://github.com/doxeylab/learn-genomics-in-linux/blob/master/task8/pheatmap.png)






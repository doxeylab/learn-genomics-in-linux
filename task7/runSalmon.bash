#!/bin/bash

for samp in SRR8451881 SRR8451882 SRR8451883 SRR8451884 SRR8451885 SRR8451886 SRR8451887 SRR8451888;
do
echo "Processing sample $samp"
salmon quant -i /fsys1/data/task4/gencode_v29_idx -l A \
         -1 /fsys1/data/task4/${samp}_1.fastq.gz \
         -2 /fsys1/data/task4/${samp}_2.fastq.gz \
         -p 6  -o quants/${samp}_quant
done

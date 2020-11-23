#!/bin/bash

for samp in SRR8451881 SRR8451882 SRR8451883 SRR8451884 SRR8451885 SRR8451886 SRR8451887 SRR8451888;
do
echo "Processing sample $samp"
salmon quant -i gencode_v29_idx -l A \
         -1 data/${samp}_1.fastq.gz \
         -2 data/${samp}_2.fastq.gz \
         -p 6  -o quants/${samp}_quant
done

#!/usr/bin/bash
RES=/public/home/lidean/NAD-RNA-seq_20211008/results
REF=/public/Reference/mouse/annotation/gencode.vM23.bed12
for id in G{1..16}; do 
    geneBody_coverage.py -r ${REF} -i ${RES}/star/R*-${id}_*align/R*-${id}_*.sorted.bam -o ${RES}/rseqc/${id}
done

#!/usr/bin/bash
TMPDIR=/public/home/lidean/NAD-RNA-seq_20211008/tmp
OUTDIR=/public/home/lidean/NAD-RNA-seq_20211008/results
ALIGNLOC=/public/home/lidean/NAD-RNA-seq_20211008/results/star
mkdir ${TMPDIR}
touch ${TMPDIR}/geneNum.txt
for bam in `ls ${ALIGNLOC}/*/*bam`; do
        id=`basename $bam .sorted.bam`                                                                                                        
        echo 'Processing ' ${bam}
        for i in {0.01,0.05,0.10,0.25,0.50,0.75,0.90}; do
                echo -e '\tSubsampling ' $i
                cat <(samtools view -H $bam) <(samtools view -@ 48 -q 255 -s $i $bam ) > ${TMPDIR}/output${i}.sam
                featureCounts -p -T 48 -a /public/Reference/mouse/annotation/gencode.vM23.annotation.gtf -o ${TMPDIR}/out${i}_counts. txt ${TMPDIR}/output${i}.sam
                geneNum=`awk '$7>10' ${TMPDIR}/out${i}_counts.txt |wc -l`
                echo -e "$i\t$geneNum\t$id" >> ${TMPDIR}/geneNum.txt
                rm -rf ${TMPDIR}/out*
        done
        TotalGeneNum=`awk '$7>10' ${OUTDIR}/featurecounts/${id}_counts/${id}_counts.txt |wc -l`
        echo -e "1.00\t$TotalGeneNum\t$id" >> ${TMPDIR}/geneNum.txt
done
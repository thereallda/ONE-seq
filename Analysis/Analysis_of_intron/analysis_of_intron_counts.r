# Analysis of intron counts 
# 13 Dec 2021
# Author: Dean Li

##-- Loading required packages --##
library(tidyverse)
library(DESeq2)
library(patchwork)
library(wesanderson)
gtf <- rtracklayer::import('gencode.vM23.annotation.gtf')

# load data from `NAD-RNA_identification.R`
load('data/NADRNA.RData')
# Intron counts ----
load("data/counts/Intron.RData")

# keep genes with intron signals (read counts > 3 in more than 50% libraries) 
intron.counts <- counts$counts
intron.counts <- intron.counts[rowSums(intron.counts > 3) >= 0.5*ncol(intron.counts),]

# Perform TPM normalization for intron ----
intron.length <- counts$annotation
intron.length <- intron.length[intron.length$GeneID %in% rownames(intron.counts),]

intron.tpm <- do.call(cbind, lapply(1:ncol(intron.counts), function(i) {
  rate = log(intron.counts[,i]) - log(intron.length$Length)
  denom = log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}))
colnames(intron.tpm) <- colnames(intron.counts)

# perform DESeq2 procedure to determine intron-enrichment ----
# intron-enrichment: fold-change >= 1 & padj < 0.05 between enrichment and input intron counts
samples_group <- c(rep('Young.Input',4),rep('Old.Input',4),rep('Young.Enrich',4),rep('Old.Enrich',4))
coldatas <- data.frame(row.names = colnames(intron.counts), group_list = samples_group)
dds.intron <- DESeqDataSetFromMatrix(countData=intron.counts, colData=coldatas, design=~group_list)
de.intron <- DESeq(dds.intron)
res.intron.y <- results(de.intron, contrast=c('group_list', 'Young.Enrich', 'Young.Input'))
res.intron.o <- results(de.intron, contrast=c('group_list', 'Old.Enrich', 'Old.Input'))

# create data.frame that holds average TPM for each group 
intron.avg <- data.frame(`Y.In`=rowMeans(intron.tpm[,1:4]),
                         `O.In`=rowMeans(intron.tpm[,5:8]),
                         `Y.En`=rowMeans(intron.tpm[,9:12]),
                         `O.En`=rowMeans(intron.tpm[,13:16]))
intron.avg$GeneID <- rownames(intron.avg)
# determine if a gene is intron-enrichment
intron.avg$Intron.stat.y <- ifelse(intron.avg$GeneID %in% rownames(subset(res.intron.y, padj < 0.05 & log2FoldChange > 0)), 'Intron.enrich', 'Intron.not.enrich')
intron.avg$Intron.stat.o <- ifelse(intron.avg$GeneID %in% rownames(subset(res.intron.o, padj < 0.05 & log2FoldChange > 0)), 'Intron.enrich', 'Intron.not.enrich')

# visualization of Intron TPM with scatter plot ----
## young ----
p1 <- ggplot(subset(intron.avg, GeneID %in% rownames(NAD_Ysig)), 
       aes(x=log2(Y.In + 1), y=log2(Y.En + 1))) +
  geom_point(aes(color=Intron.stat.y)) +
  geom_text(x=10, y=6, label='y = x') +
  ggrepel::geom_text_repel(data=subset(intron.avg, GeneID == 'ENSMUSG00000060613'), 
                           mapping=aes(label=GeneID), color='#BF351F', min.segment.length = 0) +
  labs(x='Log2(Intron TPM+1)\nInput', y='Enrichment\nLog2(Intron TPM+1)', 
       title='NAD-RNA', color='',
       subtitle=paste0('( % of Intron enrich = ', 
                       round(df1[1,]$value,3)*100, 
                       '% )')
       )
p2 <- ggplot(subset(intron.avg, !GeneID %in% rownames(NAD_Ysig)), 
             aes(x=log2(Y.In + 1), y=log2(Y.En + 1))) +
  geom_point(aes(color=Intron.stat.y)) +
  geom_text(x=12, y=5, label='y = x') +
  ggrepel::geom_text_repel(data=subset(intron.avg, GeneID == 'ENSMUSG00000089873'), 
                           mapping=aes(label=GeneID), color='#BF351F', min.segment.length = 0) +
  labs(x='Log2(Intron TPM+1)\nInput', y='Enrichment\nLog2(Intron TPM+1)',  
       title='Non-NAD-RNA', color='', 
       subtitle=paste0('( % of Intron enrich = ', 
                       round(df1[2,]$value,3)*100, 
                       '% )')
       )
ps1 <- p1 + p2 &
  geom_abline(slope=1, color='black') &
  theme_classic() &
  theme(legend.position='none',
        plot.subtitle = element_text(face='italic', color=wes_palette("Royal1", 4)[4]),
        axis.text = element_text(color='black')) &
  scale_color_manual(values=wes_palette("Royal1", 4)[c(4,1)]) 
  ggsave('results/stats/IntronEnrichment_Young.pdf', width=12, height=5)

# old ---- 
p1 <- ggplot(subset(intron.avg, GeneID %in% rownames(NAD_Osig)), 
             aes(x=log2(O.In+1), y=log2(O.En+1))) +
  geom_point(aes(color=Intron.stat.o)) +
  geom_text(x=10, y=7, label='y = x') +
  ggrepel::geom_text_repel(data=subset(intron.avg, GeneID == 'ENSMUSG00000025347'), 
                           mapping=aes(label=GeneID), color='#BF351F', min.segment.length = 0) +
  labs(x='Log2(Intron TPM+1)\nInput', y='Enrichment\nLog2(Intron TPM+1)', 
       title='NAD-RNA', color='',
       subtitle=paste0('( % of Intron enrich = ', 
                       round(df1[1,]$value,3)*100, 
                       '% )')
  ) +
  scale_y_continuous(labels=c(0,5,10,15), limits=c(0,15)) +
  scale_x_continuous(labels=c(0,5,10,15), limits=c(0,15))
p2 <- ggplot(subset(intron.avg, !GeneID %in% rownames(NAD_Ysig)), 
             aes(x=log2(O.In+1), y=log2(O.En+1))) +
  geom_point(aes(color=Intron.stat.o)) +
  geom_text(x=10, y=7, label='y = x') +
  ggrepel::geom_text_repel(data=subset(intron.avg, GeneID == 'ENSMUSG00000078673'), 
                           mapping=aes(label=GeneID), color='#BF351F', min.segment.length = 0) +
  labs(x='Log2(Intron TPM+1)\nInput', y='Enrichment\nLog2(Intron TPM+1)', 
       title='Non-NAD-RNA', color='', 
       subtitle=paste0('( % of Intron enrich = ', 
                       round(df1[2,]$value,3)*100, 
                       '% )')
  )
ps1 <- p1 + p2 &
  geom_abline(slope=1, color='black') &
  theme_classic() &
  theme(legend.position='none',
        plot.subtitle = element_text(face='italic', color=wes_palette("Royal1", 4)[4]),
        axis.text = element_text(color='black')) &
  scale_color_manual(values=wes_palette("Royal1", 4)[c(4,1)]) 

ps1
ggsave('results/stats/IntronEnrichment_Old.pdf', width=12, height=5)

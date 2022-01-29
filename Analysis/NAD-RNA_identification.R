# NAD-RNA_identification
# Merge counts files from featureCounts
##-- Loading required packages --##
library(tidyverse)
library(DESeq2)
library(ggrepel)

# In your local computer's project directory
# create the following directories
if (!dir.exists(results/qc) ) dir.create(results/qc, recursive=TRUE)
if (!dir.exists(results/EnrichedGenes) ) dir.create(results/EnrichedGenes, recursive=TRUE)

# retrieve gene annotation
gtf <- rtracklayer::import('ref/gencode.vM23.annotation.gtf')
gtf$gene_id <- gsub('\\..*$','',gtf$gene_id)
genes <- subset(gtf, type == 'gene')

# sample information 
samples_group <- factor(c(rep(c('Y.Input', 'O.Input'),each=4), rep(c('Y.Enrich', 'O.Enrich'),each=4)),
                        levels=c('Y.Input','O.Input','Y.Enrich','O.Enrich'))
samples_name <- paste(samples_group, c('R1', 'R2', 'R3', 'R4'), sep = '.')
samples_name <- samples_name %>% gsub('^Y','Young', .) %>% gsub('^O','Old', .)

# combining count files
# list the path to your count files
files_path <- list.files('data/counts/', pattern = 'G.*txt$',full.names = T)
# here use `\D` to match any character that is not a decimal digit 
# then remove them and order the path based on numreic order
files_path <- files_path[order(as.numeric(gsub("\\D", "", basename(files_path))))]
# `map_dfc()` will read in all files and combine them into a data.frame 
counts_df <- files_path %>%
  map_dfc(., ~ read.table(.x,sep = '\t', row.names = 1, skip = 2,
                          colClasses = c('character',rep('NULL',5),'integer')))
# rename columns with sample id
# `\d+` will match any digits and 
# `str_extract` extract the simple id with matched header  
colnames(counts_df) <- str_extract(basename(files_path), "G\\d+")
rownames(counts_df) <- gsub('\\..*$', '', rownames(counts_df))
write.csv(file='data/Counts.csv', x=counts_df)

# filter low-expressed genes
counts_keep <- counts_df[rowSums(counts_df > 1) >= 0.75*ncol(counts_df),]

# PCA
pca.tmp <- prcomp( t(DESeq2::vst(as.matrix(counts_keep))) )
#pca.tmp <- prcomp( t(edgeR::cpm(counts_keep, log=T)), scale=T )
rownames(pca.tmp$x) <- samples_name
# Sup.Fig.4C
# variance of each PC: summary(pca.tmp)$importance
as.data.frame(pca.tmp$x) %>% 
  mutate(group=samples_group) %>% 
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=group)) +
  geom_text_repel(label=samples_name) +
  geom_vline(xintercept=0, color='grey80', lty=2) + 
  geom_hline(yintercept=0, color='grey80', lty=2) + 
  scale_color_manual(values=wesanderson::wes_palette('Zissou1', 4)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none', axis.text = element_text(color='black')) +
  labs(x=paste0('PC1: ', round(summary(pca.tmp)$importance[2,1],3)*100, '%'), 
       y=paste0('PC2: ', round(summary(pca.tmp)$importance[2,2],3)*100, '%'))
ggsave('results/qc/PCA_PC12_vst.pdf', width=6, height=5)  

# Differential enrichment analysis
coldatas <- data.frame(row.names = colnames(counts.keep), 
                       group_list = samples_group)
dds <- DESeqDataSetFromMatrix(countData = counts_keep, 
                              colData = coldatas, 
                              design = ~ group_list)
de <- DESeq(dds)
save(de, file='data/DEobj.RData')
# retrieve DE results by contrast
NAD_Y <- results(de, contrast = c('group_list', 'Y.Enrich', 'Y.Input'))
NAD_O <- results(de, contrast = c('group_list', 'O.Enrich', 'O.Input'))
# combine DE results with gene annotation
NAD_Y <- as.data.frame(NAD_Y) %>% 
  dplyr::mutate(GeneID=rownames(.)) %>%  
  left_join(as.data.frame(genes), by=c('GeneID'='gene_id'))
NAD_O <- as.data.frame(NAD_O) %>% 
  dplyr::mutate(GeneID=rownames(.)) %>%  
  left_join(as.data.frame(genes), by=c('GeneID'='gene_id'))
# cutoff for potential NAD-RNA: log2FC >= 1 & padj < 0.05
# remove TEC genes
NAD_Ysig <- subset(NAD_Y, padj < 0.05 & log2FoldChange >= 1 & gene_type != 'TEC')
NAD_Osig <- subset(NAD_O, padj < 0.05 & log2FoldChange >= 1 & gene_type != 'TEC') 
save(NAD_Y, NAD_Ysig, NAD_O, NAD_Osig, file='data/NADRNA.RData')
# save NAD-RNA identification results in csv
NAD_Ysig %>% 
  dplyr::select(GeneID,gene_name,seqnames,gene_type,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj) %>% 
  write.csv('results/EnrichedGenes/Young_NAD-RNA.csv', row.names=F)
NAD_Osig %>% 
  dplyr::select(GeneID,gene_name,seqnames,gene_type,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj) %>% 
  write.csv('results/EnrichedGenes/Old_NAD-RNA.csv', row.names=F)

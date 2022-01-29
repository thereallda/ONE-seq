# Assessment of NudC effect
library(tidyverse)
library(DESeq2)
library(eulerr)
library(patchwork)
# read in counts file from GSM5832021 - GSM5832032
counts_df <- read.csv(file="Counts_ONE-seq_with_without_NudC.csv", row.names = 1)  

# sample information 
samples_group <- c(rep('woNudC.Enrich',2), rep('wNudC.Enrich',2), 
                   rep('woNudC.Input',2), rep('wNudC.Input',2)
                   )
samples_name <- paste(samples_group, c('1', '2'), sep = '.')

# filter low-expressed genes
counts.keep <- counts_df[rowSums(counts_df > 1) >= 0.75*ncol(counts_df),]

# Differential enrichment analysis
coldatas <- data.frame(row.names = colnames(counts.keep), 
                       group_list = samples_group)
dds <- DESeqDataSetFromMatrix(countData = counts.keep, 
                              colData = coldatas, 
                              design = ~ group_list)
de <- DESeq(dds)
save(de, file='data/DEobj.RData')

# retrieve DE results by contrast
NAD_woNudC <- results(de, contrast = c('group_list', 'woNudC.Enrich', 'woNudC.Input'), tidy=T) %>% 
  dplyr::rename(GeneID=row)
NAD_wNudC <- results(de, contrast = c('group_list', 'wNudC.Enrich', 'wNudC.Input'), tidy=T) %>% 
  dplyr::rename(GeneID=row)

# cutoff for potential NAD-RNA: log2FC >= 1 & padj < 0.05
NAD_woNudCsig <- subset(NAD_woNudC, padj < 0.05 & log2FoldChange >= 1)
NAD_wNudCsig <- subset(NAD_wNudC, padj < 0.05 & log2FoldChange >= 1) 
save(NAD_woNudC, NAD_woNudCsig, NAD_wNudC, NAD_wNudCsig, file='data/NADRNA.RData')

# Merge DE results
NAD_all <- full_join(NAD_woNudC, NAD_wNudC, by='GeneID', suffix=c('.wo', '.w'))

# All identified NAD-RNA  
NAD_ID <- union(NAD_wNudCsig$GeneID, NAD_woNudCsig$GeneID)

# NudC- specific genes
woNudC_ID <- subset(NAD_woNudCsig, !GeneID %in% NAD_wNudCsig$GeneID)$GeneID

# Visualization
# venn diagram
vp <- plot(euler(list(
  'NudC-' = NAD_woNudCsig$GeneID,
  'NudC+' = NAD_wNudCsig$GeneID
  )), 
  edges=NA, 
  labels=c(paste0('NudC-\n',length(woNudC_ID)),
           paste0('NudC+\n',length(NAD_ID)-length(woNudC_ID))),
  fill=c('#666666', '#C8B98C', '#C8B98C'))

# scatter plot
sp <- NAD_all %>% 
  filter(GeneID %in% NAD_ID) %>% 
  mutate(woNudC.specific = if_else(GeneID %in% woNudC_ID, 'Expelled', 'Within')) %>% 
  ggplot(aes(x = log2FoldChange.w,
             y = log2FoldChange.wo,
             color = woNudC.specific)) +
  geom_point() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'),
        legend.position = 'none') +
  scale_color_manual(values=c('#666666', '#C8B98C')) +
  labs(x='Log2 Fold Change (Enrichment/Input)\nNudC+', 
       y='NudC-\nLog2 Fold Change (Enrichment/Input)')

wrap_elements(vp)/sp + plot_layout(heights = c(1, 2))
ggsave('results/NudC_assessment.pdf', width=6, height=8)

# Assessment of spike-in RNA
library(tidyverse)
library(patchwork)

# read in counts file from GSM5831997 - GSM5832020
counts_df <- read.csv(file="Counts_ONE-seq_Spikein_RNA.csv", row.names = 1)  

# sample information 
samples_group <- c(rep('NAD0.Input',3), rep('NAD0.Enrich',3), 
                   rep('NAD1.Input',3), rep('NAD1.Enrich',3),
                   rep('NAD5.Input',3), rep('NAD5.Enrich',3),
                   rep('NAD10.Input',3), rep('NAD10.Enrich',3)
                   )
samples_name <- paste(samples_group, c('1', '2', '3'), sep = '.')

# filter low-expressed genes
counts.keep <- counts_df[rowSums(counts_df>0) > 0.75*ncol(counts_df), ]
# perform CPM normalization followed by log transformation 
lcpm.keep <- edgeR::cpm(counts.keep, log=T)

# ComBat for removing batch effect
# input clean and pre-normalized data
lcpm.rmbat <- sva::ComBat(lcpm.keep, batch=rep(c(1,2), each=12))

# Calculate the logFoldChange of the SpikeinRNA CPM between Enrichment and Input
# NAD 0%, 1%, 5%, 10%
lfc <- c(lcpm.rmbat['SpikeinRNA', 4:6]-lcpm.rmbat['SpikeinRNA', 1:3], 
         lcpm.rmbat['SpikeinRNA', 10:12]-lcpm.rmbat['SpikeinRNA', 7:9],
         lcpm.rmbat['SpikeinRNA', 16:18]-lcpm.rmbat['SpikeinRNA', 13:15],
         lcpm.rmbat['SpikeinRNA', 22:24]-lcpm.rmbat['SpikeinRNA', 19:21])
lfc <- setNames(object = lfc, nm=rep(paste0('NAD', c(0,1,5,10)), each=3))

# boxplot visualization
tibble(
  LFC=lfc,
  NAD=factor(names(lfc), levels=c('NAD0','NAD1','NAD5','NAD10'))
) %>% 
  ggplot(aes(x=NAD, y=2^LFC)) +
  geom_boxplot(width=0.5) +
  geom_point() +
  geom_hline(yintercept = 2, color='#C93312', lty='dashed') +
  theme_minimal() +
  theme(axis.ticks.y = element_line(color='black'),
        axis.text = element_text(color='black'),
        axis.line.y.left = element_line(color='black'),
        axis.line.x.bottom = element_line(color='black'),
        panel.grid = element_blank()
  ) +
  scale_x_discrete(labels=c('0%', '1%', '5%', '10%')) +
  labs(x='Ratio of NAD-RNA', 
       y='Spike-in\nFold change (Enrichment/Input)')
ggsave('results/spike-in_FC.pdf', width=4, height=3)
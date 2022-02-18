# NAD-RNA_pathways
library(gprofiler2)
library(tidyverse)
load('data/NADRNA.RData')

# In your local computer's project directory
# create the following directories
if (!dir.exists(results/FE) ) dir.create(results/FE, recursive=TRUE)

# sources: GO:BP and REAC only
multi_gores <- gost(query=list('Young'=NAD_Ysig$GeneID, 'Old'=NAD_Osig$GeneID),
                     organism='mmusculus',
                     exclude_iea=TRUE, 
                     correction_method = 'gSCS', # default 'gSCS'
                     evcodes=TRUE,
                     sources=c('GO:BP', 'REAC'))
# filtering pathways with less than 5 genes and more than 350 genes and 
# sharing gene with queries less than 3 genes
multi_gores_cut <- list("result" = subset(multi_gores$result, 
                                           term_size > 5 & 
                                             term_size < 350 & 
                                             intersection_size > 3),
                         "meta" = multi_gores$meta)
# convert ENSEMBL id to SYMBOL
multi_gores_cut$result$intersection %<>% 
  str_split(',') %>% 
  map_chr(function(x) subset(NAD_Ysig, GeneID %in% x)$gene_name %>% paste(collapse=','))

# Creating a Generic Enrichment Map (GEM) file
gem <- multi_gores2_cut$result[,c("term_id", "term_name", "p_value", "intersection", "query")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes", "Phenotype")
gem$FDR <- gem$p.Val
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
# output young and old gem files for visualization of pathway network in cytoscape, respectively
subset(gem, Phenotype == 'Young') %>% 
  write.table(file = "results/FE/YO_GOBP_REAC_gemY.txt", sep="\t", quote=F, row.names=F)
  
subset(gem, Phenotype == 'Old') %>% 
  write.table(file = "results/FE/YO_GOBP_REAC_gemO.txt", sep="\t", quote=F, row.names=F)

# bar plot of top 10 go terms in young
multi_gores_cut$result %>% 
  filter(query == 'Young') %>%
  filter(term_name != 'Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.') %>%
  slice_min(p_value,n=10) %>%
  select(term_name, intersection_size, p_value) %>% 
  mutate(p_value = -log10(p_value),
         term_name = str_to_sentence(term_name)) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = factor(term_name, levels=rev(unique(term_name))))) +
  geom_bar(aes(y=value, fill=variable), stat='identity', 
           position=position_dodge(width=0.9), width=0.8) +
  coord_flip() +
  theme_minimal() +
  geom_hline(yintercept = (-log10(0.05)),lty=4,col="grey",lwd=0.6) +
  theme(axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(color='black'),
        axis.line = element_line(),
        axis.ticks = element_line(color='black'),
        legend.position = 'top',
        panel.grid = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_manual(values = c("#8696A1","#CDDEE9"), 
                    labels = c('Gene Counts', '-log10(p.adjust)')) +
  labs(x='',y='',fill='')
ggsave('results/FE/Pathway_Young_top10.pdf', width=8, height=4)

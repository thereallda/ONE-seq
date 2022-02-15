# NAD-RNA_characterization
##-- Loading required packages --##
library(tidyverse)
library(patchwork)
library(ggside)
library(ggrepel)
library(gghalves)
library(ggstatsplot)
library(wesanderson)
# load data from `NAD-RNA_identification.R`
load('data/DEobj.RData')
load('data/NADRNA.RData')

# In your local computer's project directory
# create the following directories
if (!dir.exists(results/stats) ) dir.create(results/stats, recursive=TRUE)

# make object for visualization with ggplot2
nadtmp1 <- tibble(
  GeneID=c(NAD_Ysig$GeneID, NAD_Osig$GeneID),
  LFCs=c(NAD_Ysig$log2FoldChange, NAD_Osig$log2FoldChange),
  Groups=factor(c(rep('Young',nrow(NAD_Ysig)), rep('Old',nrow(NAD_Osig))), 
                  levels=c('Young','Old')
                )
  )
genes_anno <- NAD_Y[,c('GeneID','gene_name','seqnames','gene_type')]
ggdf1 <- inner_join(nadtmp1, genes_anno, by='GeneID')
# rename all *pseudogene as pseudogene
# consider lncRNA, snRNA, miRNA, scaRNA as ncRNA
# consider IG_C_gene, TR_C_gene as protein_coding
ggdf1$gene_type <- str_replace(ggdf1$gene_type, '.*pseudogene', 'pseudogene')
ggdf1$gene_type <- str_replace(ggdf1$gene_type, 'lncRNA|snRNA|miRNA|scaRNA', 'ncRNA')
ggdf1$gene_type <- str_replace(ggdf1$gene_type, 'IG_C_gene|TR_C_gene', 'protein_coding')
ggdf1$gene_type <- factor(ggdf1$gene_type, levels=c('ncRNA','pseudogene','protein_coding'))

# MA plot ----
# young
NAD_Y %>% 
  mutate(sig=if_else(log2FoldChange > 1 & padj < 0.05, 'NAD-RNA', 'None')) %>% 
  ggplot(aes(x=log2(baseMean+1), y=log2FoldChange)) +
  geom_point(aes(color=sig), alpha=0.4, shape=19) +
  annotate('text', x=19, y=3, label='bold("2017")', color="#046C9A", size=5, parse=T) +
  theme_classic() +
  theme(axis.text = element_text(color='black'),
        axis.ticks = element_line(color='black'),
        legend.position = 'top') +
  geom_hline(yintercept = 0, lty='dashed') +
  scale_color_manual(values = c("#046C9A", "#D5D5D3")) +
  labs(x='Log2 (average abundance+1)', y='Log2 Enrichment/Input', color='', title='Young')
ggsave('results/stats/MAplot_young.pdf', width=5, height=4)

# old
NAD_O %>% 
  mutate(sig=if_else(log2FoldChange > 1 & padj < 0.05, 'NAD-RNA', 'None')) %>% 
  ggplot(aes(x=log2(baseMean+1), y=log2FoldChange)) +
  geom_point(aes(color=sig)) +
  annotate('text', x=19, y=3, label='bold("1820")', color="#046C9A", size=5, parse=T) +
  theme_classic() +
  theme(axis.text = element_text(color='black'),
        axis.ticks = element_line(color='black'),
        legend.position = 'top') +
  geom_hline(yintercept = 0, lty='dashed') +
  scale_color_manual(values = c("#046C9A", "#D5D5D3")) +
  labs(x='Log2 (average abundance+1)', y='Log2 Enrichment/Input', color='', title='Old')
ggsave('results/stats/MAplot_old.pdf', width=5, height=4)

# gene type ----
ggdf1 %>% 
  dplyr::count(Groups, gene_type) %>% 
  group_by(Groups) %>% 
  mutate(pct=n/sum(n)) %>% 
  ggplot(aes(x = Groups, y = n, group=gene_type)) +
  geom_bar(aes(fill=gene_type), stat='identity', width = 0.5) +
  geom_label(aes(label=paste0(n,"\n(",sprintf("%1.1f", pct*100),"%)")), 
             position=position_stack(vjust=0.5), size=3) +
  scale_fill_manual(values = c("#ECCBAE", "#D69C4E", "#046C9A")) +
  scale_x_discrete(labels=c(paste0('Young\n(n=',sum(ggdf1$Groups == 'Young'),')'),
                            paste0('Old\n(n=',sum(ggdf1$Groups == 'Old'),')')
  )) +
  theme_classic() +
  theme(axis.text = element_text(color='black', size=12),
        legend.position = 'top') +
  labs(x = '', y = '',fill = '')
ggsave('results/stats/NAD-RNA_Type.pdf', width=4, height=6)

# chromosome distribution ----
hp1 <- as.data.frame(table(NAD_Ysig$seqnames)/table(NAD_Y$seqnames)*100) %>% 
  ggplot(aes(x=Var1, y=round(Freq, 2))) +
  geom_bar(stat='identity', fill="#798E87", width=0.8) +
  geom_text(aes(label=round(Freq, 1), y=1.05*round(Freq, 1)), size=3) +
  theme_classic() +
  theme(axis.text = element_text(color='black'),
        axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.5),
        axis.ticks.y = element_line(color='black'),
        axis.ticks.x = element_blank()) + 
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x='', y='Proportion of NAD-RNAs', title='Young')

hp2 <- as.data.frame(table(NAD_Osig$seqnames)/table(NAD_O$seqnames)*100) %>% 
  ggplot(aes(x=Var1, y=round(Freq, 2))) +
  geom_bar(stat='identity', fill="#798E87", width=0.8) +
  geom_text(aes(label=round(Freq, 1), y=1.05*round(Freq, 1)), size=3) +
  theme_classic() +
  theme(axis.text = element_text(color='black'),
        axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.5),
        axis.ticks.y = element_line(color='black'),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x='', y='Proportion of NAD-RNAs', title='Old')
hp1 + hp2
ggsave('results/stats/ChrDis_NAD-RNA.pdf', width=14, height=4)

# gene-length by fold change group ----
# retrieve gene length from featrueCounts files
Translen <- read.table('data/S1_counts.txt', sep='\t', header=T)
ggdf2 <- inner_join(nadtmp1, Translen[,c('Geneid','Length')], by = c('GeneID'='Geneid')) %>% 
  group_by(Groups) %>% 
  mutate(quartile=ntile(LFCs, 10)) %>%
  ungroup()
# half box half point
ggplot(ggdf2, aes(x=factor(quartile), y=Length/1000)) +
  geom_half_boxplot(aes(fill=factor(quartile)), color="grey30",
                    outlier.shape = NA) +
  geom_half_point_panel(aes(color=LFCs), side = 'r',size=1.5,alpha=.6, stroke=0,
                        position=position_jitter(width=.15)) +
  see::theme_modern() +
  theme(axis.ticks = element_line(color='black'),
        axis.text = element_text(color='black'),
        legend.position='right') +
  facet_wrap(.~Groups, nrow=2, scales = 'free_y') +
  scale_fill_manual(values=colorRampPalette(c('#e5e5be','#C8B98C',"#567794",'#38678F'))(10), guide='none') +
  scale_color_gradientn(colors=colorRampPalette(c('#e5e5be','#C8B98C',"#567794",'#38678F'))(30),
                       breaks=1:3, labels=1:3, limits=c(0,3.5)) +
  labs(x='Deciles', y='Gene Length(kb)',color='log2FC')
ggsave('results/stats/FC_GeneLength.pdf', width=8, height=6)

# Dynamics of global NAD modification levels in violin-box-scatter plot ----
ggbetweenstats(ggdf1, x=Groups, y=LFCs,
               ggtheme = theme_classic()) +
  theme(axis.text = element_text(color='black')) +
  scale_color_manual(values=c('#38678F', '#C8B98C')) +
  scale_y_continuous(labels=1:4, breaks=1:4, limits=c(0.8,4),expand=c(0,0)) +
  labs(x='', y='log2 Fold Change (Enrichment/Input)')
ggsave('results/stats/NAD-RNA_Ratio.pdf', width=5, height=6)

# Dynamics of gene-specific NAD modification levels in scatter plot ----
cols1 <- c('GeneID','gene_name','log2FoldChange')
nadID <- union(NAD_Ysig$GeneID, NAD_Osig$GeneID)

NAD_YO <- inner_join(NAD_Y[,cols1], NAD_O[,cols1], by='GeneID') %>% 
  filter(GeneID %in% nadID) %>% 
  mutate(FCDiff_OY=2^log2FoldChange.y/2^log2FoldChange.x) %>% 
  mutate(Condition=ifelse(FCDiff_OY > 3/2, 
                          yes='Old', 
                          no=ifelse(FCDiff_OY < 2/3, 
                                    yes='Young', no='Common')),
         Condition=factor(Condition, levels=c('Young','Common','Old')))

ggplot(NAD_YO, aes(x=2^log2FoldChange.x, y=2^log2FoldChange.y)) +
  geom_point(aes(color=Condition)) +
  geom_segment(aes(x=0,xend=2, y=2, yend=2), linetype='dashed') +
  geom_segment(aes(x=2,xend=2, y=0, yend=2), linetype='dashed') +
  geom_abline(slope = 2/3, linetype='dashed', color='grey20', size=1) +
  geom_abline(slope = 3/2, linetype='dashed', color='grey20', size=1) +
  geom_abline(slope = 1, linetype='dashed', color='#C93312', size=1) +
  scale_color_manual(values = c('#38678F', 'grey', '#C8B98C')) +
  annotate('text', x=9.5, y=0.5, label="247", color='#38678F') +
  annotate('text', x=1, y=8.5, label="105", color='#C8B98C') +
  annotate('text', x=9, y=8, label='1226', color='grey') +
  annotate('text', x=6.8, y=8.5, label='772', color='grey') +
  geom_xsidedensity(aes(y=stat(density)),fill='#52616a') +
  geom_ysidedensity(aes(x=stat(density)),fill='#52616a') +
  theme_bw() + 
  theme(axis.text=element_text(color='black'),
        legend.position='bottom',
        panel.grid=element_blank(), 
        ggside.panel.border = element_blank(),
        ggside.axis.text = element_blank(),
        ggside.axis.ticks = element_blank()
        ) +
  scale_x_continuous(labels = 0:10, breaks = 0:10, limits = c(0,10), expand = c(0,0)) +
  scale_y_continuous(labels = 0:10, breaks = 0:10, limits = c(0,10), expand = c(0,0)) +
  scale_xsidey_continuous(limits = c(0,1)) +
  scale_ysidex_continuous(limits = c(0,1)) +
  labs(color='', x='Fold Change (Enrichment/Input)\nYoung', y='Old\nFold Change (Enrichment/Input)', title='All NAD-RNA')
ggsave('results/stats/YO_FCplot.pdf', width=6, height=6)
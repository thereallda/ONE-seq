# seuqncing stats
library(tidyverse)
# sequencing saturation
samples_group <- factor(c(rep(c('Young.Input', 'Old.Input'),each=4), 
                          rep(c('Young.Enrich', 'Old.Enrich'),each=4)),
                        levels=c('Young.Input','Old.Input','Young.Enrich','Old.Enrich'))
geneNums <- read.table('results/qc/geneNum.txt', sep='\t',header=F)
colnames(geneNums) <- c('Per','Num','Sample')

sample.labs <- set_names(x=paste(samples_group, paste0('R',1:4), sep='.'), 
                         nm=paste0('G',1:16))
geneNums %>% 
  mutate(Sample=factor(str_extract(Sample, "G\\d+"), levels=paste0('G',1:16))) %>% 
  add_row(Per=rep(0,16), 
          Num=rep(0,16), 
          Sample=factor(paste0('G',1:16), levels=paste0('G',1:16))) %>% 
  arrange(Sample, Per) %>% 
  ggplot(aes(x=Per, y=Num)) +
  geom_line() +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line('black'),
        # strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold'),
        axis.text = element_text(color='black')) +
  facet_wrap(~Sample, labeller = labeller(Sample=sample.labs)) +
  scale_x_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(x='Fraction of library', y = 'Numbers of Mapped Genes')
ggsave('results/qc/SeqSaturation.pdf', width = 7.5, heigh = 6)
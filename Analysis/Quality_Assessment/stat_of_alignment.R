# alignment stat
library(tidyverse)
library(ggh4x)
# unique mapped read pairs
maprate <- read.table('results/qc/star_alignment.tsv', sep='\t',header=T,check.names=F)
colnames(maprate) <- c('Category','uni','multi1','multi2','unmap1','unmap2')

maprate.sorted <- maprate %>% 
  mutate(ID=factor(str_extract(Category, "G\\d+"), levels=paste0('G',1:16))) %>% 
  arrange(ID) %>% 
  mutate(age=str_split(samples_group, '\\.', simplify=T)[,1],
         assay=str_split(samples_group, '\\.', simplify=T)[,2],
         group=paste0(age,c('_R1.','_R2.','_R3.','_R4.'), assay)
         ) 

ggplot(maprate.sorted, aes(x=weave_factors(group), y=uni/1e6)) +
  geom_bar(stat="identity", width=0.8, fill='grey60') +
  geom_text(aes(label=round(uni/1e6,1)), hjust=-0.2) +
  coord_flip() + 
  theme_bw() +
  scale_x_discrete(guide='axis_nested', 
                   limits=rev(levels(weave_factors(maprate.sorted$group))),
                   labels=rev(gsub('_',' ',maprate.sorted$group))) +
  theme(axis.text = element_text(color='black'),
        ggh4x.axis.nestline.y = element_line(size=1),
        ggh4x.axis.nesttext.y = element_text(color='black', face='bold', size=12),
        panel.grid = element_blank()) +
  labs(x='',y='Qualified read pairs (M)')
ggsave('results/qc/UniqueMap_lib.pdf', width=8, height=4)

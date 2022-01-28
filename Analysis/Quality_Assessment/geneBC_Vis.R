# gene body coverage visualization in R
library(tidyverse)
# read in all coverage data and merage as one data.frame
files_path <- list.files('results/GeneBodyCoverage/', pattern='txt$', full.names=T)
files_path <- files_path[order(as.numeric(gsub("\\D", "", basename(files_path))))]
coverage_df <- files_path %>% 
  map_dfc(function(i) data.table::transpose(read.table(i, sep='\t', header=T), make.names=1))
colnames(coverage_df) <- str_extract(basename(files_path), "G\\d+")

tmp1 <- t( t(coverage_df) - apply(coverage_df, 2, min) )
tmp2 <- t( t(tmp1)/apply(tmp1, 2, max) )
coverage_df <- as.data.frame(tmp2)

samples_group <- factor(c(rep(c('Y.Input', 'O.Input'),each=4), rep(c('Y.Enrich', 'O.Enrich'),each=4)),
                        levels=c('Y.Input','O.Input','Y.Enrich','O.Enrich'))
samples_df <- data.frame(ID=paste0('G',1:16), Groups=samples_group)

ggdf1 <- coverage_df %>%
  reshape2::melt() %>% 
  mutate(Percentile=rep(1:100, times=ncol(coverage_df))) %>% 
  left_join(samples_df, by=c('variable'='ID'))
colnames(ggdf1) <- c('Samples', 'Coverage', 'Percentile', 'Groups')
ggdf1$Age <- factor(str_split(ggdf1$Groups, '\\.', simplify = T)[,1], levels=c('Y','O'))
ggdf1$Assay <- factor(str_split(ggdf1$Groups, '\\.', simplify = T)[,2], levels=c('Input','Enrich'))

sample_label <- c(
  `Y` = 'Young',
  `O` = 'Old'
)

ggplot(ggdf1, aes(x=Percentile, y=Coverage, group=Samples, color=Groups)) +
  geom_line() +
  facet_wrap(~Age, labeller=as_labeller(sample_label)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(color='black'),
        strip.background = element_rect(fill=NA, linetype='solid')) +
  scale_color_manual(values=wesanderson::wes_palette('Zissou1')) +
  labs(x="Gene body percentile (5'->3')",color='')
ggsave('results/GeneBodyCoverage/GBC_Vis_v1.pdf', height = 4, width = 8)

# New facet label names for Age variable
age.labs <- c("Young", "Old")
names(age.labs) <- c("Y", "O")

# New facet label names for Assay variable
assay.labs <- c("Input", "Enrich")
names(assay.labs) <- c("Input", "Enrich")

ggplot(ggdf1, aes(x=Percentile, y=Coverage, group=Samples, color=Groups)) +
  geom_line() +
  facet_grid(Age~Assay, labeller = labeller(Age=age.labs, Assay=assay.labs)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(color='black'),
        strip.background = element_rect(fill=NA, linetype='solid')) +
  scale_color_manual(values=wesanderson::wes_palette('Zissou1')) +
  labs(x="Gene body percentile (5'->3')",color='')
ggsave('results/GeneBodyCoverage/GBC_Vis_v2.pdf', height = 6, width = 8)

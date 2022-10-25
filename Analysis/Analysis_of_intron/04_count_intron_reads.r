# Counting intron reads
####
# Source: https://github.com/charitylaw/Intron-reads/blob/master/Analyses/analysis_of_intron_exploration.R
# 13 Dec 2021 (modified by Dean Li)
# 10 Jan 2017 (Last updated 7 Jun 2018)
# Charity Law and Albert Zhang

# Fields to define
exon.anno <- "data/ref/Exon.txt" # an exon saf file
genebody.anno <- "data/ref/Genebody.txt" # a genebody saf file
bam <- list.files('data/bam', full.names = TRUE) # bam file names
bam <- bam[order(as.numeric(gsub("\\D", "", basename(bam))))] # order bam by sample id 
is.paired.end <- TRUE # seq protocol for read summarisation
sample.names <- paste0('G', 1:16) # sample names (unique)
groups <- sample.names # groups for samples (non-unique)

# Load libraries
library(limma)
library(edgeR)
library(Rsubread)
library(stringr)

# Create folders
if (!dir.exists('data/counts')) dir.create("data/counts", showWarnings=FALSE)
if (!dir.exists('results')) dir.create("results", showWarnings=FALSE)
if (!dir.exists('results/plots')) dir.create("results/plots", showWarnings=FALSE)


### SUMMARISE READS TO COUNTS

# Exon
anno <- read.delim(exon.anno)
# order exon annotation by gene id 
anno <- anno[order(anno$GeneID),]
counts <- featureCounts(annot.ext=anno, files=bam, isPairedEnd=is.paired.end, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, strandSpecific=0, nthreads=16)
colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
save(counts, file="data/counts/Exon.RData")  
# Genebody
anno <- read.delim(genebody.anno)
counts <- featureCounts(annot.ext=anno, files=bam, isPairedEnd=is.paired.end, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, strandSpecific=0, nthreads=8)
colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
save(counts, file="data/counts/Genebody.RData")    
# Intron  
load("data/counts/Exon.RData")
exon <- counts
load("data/counts/Genebody.RData")
genebody <- counts
counts$counts <- pmax(genebody$counts-exon$counts,0)
save(counts, file="data/counts/Intron.RData")

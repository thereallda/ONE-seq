# Framework of NAD-RNA sequencing data analysis



## RNA-seq_pipeline

In "RNA-seq_pipeline" folder:

- `STAR_pipeline.sh`: In-house script for RNA-sequencing data quality control, mapping and quantification. 
- `README.md`: Detailed description of the `STAR_pipeline.sh`



## Assessment of data quality

Scripts in "Quality_Assessment" folder: 

### Quality of alignment

`star_alignment.tsv`: Statistics of star alignment.

`stat_of_alignment.R`: Summarizing and visualizing the uniquely-mapped read numbers.



### Sequencing saturation

`seqSaturation.sh`: In-house script for sequencing saturation assessment.

`seqSaturation_Vis.R`: Visualization of sequencing saturation in R.



### Assessment of RNA integrity 

`geneBC.sh`: Calculate the RNA-seq reads coverage over gene body using `RseQC`.

`geneBC_Vis.R`: Visualization of reads coverage over gene body in R.



## NAD-RNA analysis

The framework of main data analysis of NAD-RNA is suitable for any NAD-RNA-seq data processed with the RNA-seq pipeline mentioned above. Below describing the analysis with data from young and old mouse.

> Table of sample information.
>
> | GSM_ID                                                       | Sample_Name       |
> | ------------------------------------------------------------ | ----------------- |
> | [GSM5832033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832033) | Young_Input_rep1  |
> | [GSM5832034](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832034) | Young_Input_rep2  |
> | [GSM5832035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832035) | Young_Input_rep3  |
> | [GSM5832036](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832036) | Young_Input_rep4  |
> | [GSM5832037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832037) | Old_Input_rep1    |
> | [GSM5832038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832038) | Old_Input_rep2    |
> | [GSM5832039](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832039) | Old_Input_rep3    |
> | [GSM5832040](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832040) | Old_Input_rep4    |
> | [GSM5832041](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832041) | Young_Enrich_rep1 |
> | [GSM5832042](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832042) | Young_Enrich_rep2 |
> | [GSM5832043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832043) | Young_Enrich_rep3 |
> | [GSM5832044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832044) | Young_Enrich_rep4 |
> | [GSM5832045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832045) | Old_Enrich_rep1   |
> | [GSM5832046](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832046) | Old_Enrich_rep2   |
> | [GSM5832047](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832047) | Old_Enrich_rep3   |
> | [GSM5832048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832048) | Old_Enrich_rep4   |

### Identification

`NAD-RNA_identification.R`: Script for identification of NAD-RNA.

The analysis requires:

- Counts files from `featureCounts`.
- Sample information (sample names, group, contrast).
- Annotation files ([gencode.vM23.annotation.gtf](https://www.gencodegenes.org/mouse/release_M23.html)).

The analysis includes: 

- Data cleaning.
- Creating PCA plots.
- Differential enrichment analysis for NAD-RNA identification.

### Characterization

`NAD-RNA_characterization.R`: Script for exploring the characteristics of NAD-RNA.

The analysis requires:

- Results of differential enrichment analysis from `NAD-RNA_identification.R`.
- Sample information (sample names, group, contrast).
- Annotation files ([gencode.vM23.annotation.gtf](https://www.gencodegenes.org/mouse/release_M23.html)).

The analysis explores the characteristics of NAD-RNA, including: 

- gene abundance ~ NAD modification level (MA-plot);
- Gene type;
- Chromosome distribution;
- Gene length;
- Dynamics of global NAD modification levels;
- Dynamics of gene-specific NAD modification levels.



### Pathway enrichment analysis

`NAD-RNA_pathways.R`: Script for exploring the pathways associated with NAD-RNA.

The analysis requires:

- Results of differential enrichment analysis from `NAD-RNA_identification.R`.
- Sample information (sample names, group, contrast).
- `GO:BP` and `REAC` pathways from `gprofiler` (https://biit.cs.ut.ee/gprofiler/static/gprofiler_mmusculus.name.zip) 

The analysis includes:

- Pathway enrichment analysis of NAD-RNA.
- Visualization of pathway enrichment results.



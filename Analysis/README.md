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

The framework of main data analysis of NAD-RNA is suitable for any NAD-RNA-seq data processed with the RNA-seq pipeline mentioned above. Below demonstrating the analysis with data from young and old mouse.

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



## Others

Others folder contains the following scripts: 

### spike-in RNA assessment

This analysis is performed with NAD-RNA-seq data from samples mixed with spike-in RNA that had different ratio of NAD-RNA.


> Table of sample information.
> | GSM_ID                                                       | Sample_Name       |
> | ------------------------------------------------------------ | --------------- |
> | [GSM5831997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5831997) | NAD0_Input_rep1 |
> | [GSM5831998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5831998) | NAD0_Input_rep2 |
> | [GSM5831999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5831999) | NAD0_Input_rep3 |
> | [GSM5832000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832000) | NAD1_Input_rep1   |
> | [GSM5832001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832001) | NAD1_Input_rep2   |
> | [GSM5832002](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832002) | NAD1_Input_rep3   |
> | [GSM5832003](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832003) | NAD5_Input_rep1   |
> | [GSM5832004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832004) | NAD5_Input_rep2   |
> | [GSM5832005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832005) | NAD5_Input_rep3   |
> | [GSM5832006](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832006) | NAD10_Input_rep1  |
> | [GSM5832007](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832007) | NAD10_Input_rep2  |
> | [GSM5832008](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832008) | NAD10_Input_rep3  |
> | [GSM5832009](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832009) | NAD0_Enrich_rep1  |
> | [GSM5832010](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832010) | NAD0_Enrich_rep2  |
> | [GSM5832011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832011) | NAD0_Enrich_rep3  |
> | [GSM5832012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832012) | NAD1_Enrich_rep1  |
> | [GSM5832013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832013) | NAD1_Enrich_rep2  |
> | [GSM5832014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832014) | NAD1_Enrich_rep3  |
> | [GSM5832015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832015) | NAD5_Enrich_rep1  |
> | [GSM5832016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832016) | NAD5_Enrich_rep2  |
> | [GSM5832017](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832017) | NAD5_Enrich_rep3  |
> | [GSM5832018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832018) | NAD10_Enrich_rep1 |
> | [GSM5832019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832019) | NAD10_Enrich_rep2 |
> | [GSM5832020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832020) | NAD10_Enrich_rep3 |

Following analysis described above (i.e., `RNA-seq_pipeline`), the read counts of spike-in RNA were assessed. Raw counts of spike-in RNA can be accessed in the counts file online ([GSE194271_Counts_ONE-seq_Spikein_RNA.csv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194271&format=file&file=GSE194271%5FCounts%5FONE%2Dseq%5FSpikein%5FRNA%2Ecsv%2Egz)).



`spike-in_RNA_assessment.R`: Data analysis script for assessment of spike-in RNA.

The analysis requires:

- Counts files from `featureCounts`.
- Sample information (sample names, group, contrast).

The analysis includes: 

- Data cleaning.
- Creating boxplot of the normalized read counts of spike-in RNA.



### NudC assessment

This analysis is performed with NAD-RNA-seq data with or without NudC treatment.

> Table of sample information.
>
> | GSM_ID                                                       | Sample_Name              |
> | ------------------------------------------------------------ | ------------------------ |
> | [GSM5832021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832021) | Without_NudC_Input_rep1  |
> | [GSM5832022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832022) | Without_NudC_Input_rep2  |
> | [GSM5832024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832024) | With_NudC_Input_rep1     |
> | [GSM5832025](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832025) | With_NudC_Input_rep2     |
> | [GSM5832027](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832027) | Without_NudC_Enrich_rep1 |
> | [GSM5832028](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832028) | Without_NudC_Enrich_rep2 |
> | [GSM5832030](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832030) | With_NudC_Enrich_rep1    |
> | [GSM5832031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5832031) | With_NudC_Enrich_rep2    |

Following analysis described above (i.e., `RNA-seq_pipeline` and `NAD-RNA_identification`), the NAD modification level of NAD-RNA identified from the NAD-RNA-seq data with or without NudC treatment were assessed. Information of NAD-RNA can be found in our online supplementary information (Supplementary Table 1). 



`NudC_assessment.R`: Script for assessing the effect of NudC treatment.

The analysis requires:

- Counts files from `featureCounts`.
- Sample information (sample names, group, contrast).
- Results of differential enrichment analysis from `NAD-RNA_identification.R`.

The analysis includes: 

- Overlaps of NAD-RNA identified from different condition.
- Visualization with venn diagram and scatter plot.




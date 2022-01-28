# Description

The `STAR_pipeline.sh`  integrates the common steps for RNA-seq analysis, including : 

1. Trimming and quality control (`Trim Galore` and `FastQC`)

2. Mapping (`STAR`)

3. sam2bam, sort and index (`samtools`)

4. gene expression quantification (`featureCounts`)




# Dependences

```shell
python #v3.7
fastqc #v0.11.9
trim_galore #v0.6.6
STAR #2.7.6a
samtools #1.11
featureCounts #2.0.1
```



# Prerequisite 

In order to perform analysis with `STAR_pipeline.sh`, a project directory and the sub-directories must be created as the following structure: 

> Project_directory/
> |-- data
> |   -- fastq
> |   	`*fastq.gz`
> |-- results
> |-- src

```bash
mkdir -p data/fastq
mkdir results
mkdir src
```

- `data/fastq/` for holding all `fastq` files to be analysis. Put all your `fastq` files in `<project_dir>/<data_dir>/fastq/`. 

- `results/` for holding results.

- `src/` for any source files, e.g., `STAR_pipeline.sh`..



# Usage
```Shell
bash STAR_pipeline.sh [--pair] -d <project_dir> -o <output_dir> -i <input_data_dir> --ref <reference_genome> --gtf <GTF_file> -t <threads>
```

`-d`: project directory, e.g.,  `~/RNA-seq_pipe`, absolute path is preferable. 

`-i`: data directory, un-processed `fastq` files are in `<data_dir>/fastq` and trimmed `fastq` files will be generated at `<data_dir>/clean`.

`-o`: output directory, create results in the specified output directory. 

`--ref`: path to the directory where genome files are stored. 

`--gtf`:  path to the annotation files (in GTF/GFF format).  Gzipped file is also accepted.

`--pair`: `FLAG` turn on paired-end processing mode. 

`-t`: `INT`, number of threads



## Example 

First, create a project directory, for example `RNApipeTest` 

```bash
mkdir RNApipeTest
```

In `<project_dir>`, create the directories mentioned above and copy all the `fastq` files in `<data_dir>/fastq` . Also, place `STAR_pipeline.sh` at `src/`. 

```bash
cd RNApipeTest
mkdir -p data/fastq
mkdir results
mkdir src
cp /wherever/you/store/the/fastq.gz ./data/fastq
```

After preparation of the project directory, it should look like:

> RNApipeTest/
> ├── data
> │   └── fastq
> │       ├── S1_R1.fastq.gz
> │       ├── S1_R2.fastq.gz
> │       ├── S2_R1.fastq.gz
> │       ├── S2_R2.fastq.gz
> │       ├── ...
> │       └── S20_R2.fastq.gz
> ├── results
> └── src
>     └── STAR_pipeline.sh



**Before entering the following command, make sure the current working directory contains `STAR_pipeline.sh`.** In addition, genome reference index for `STAR` alignment and annotation for `featureCounts` should be prepared before running the pipeline.

```shell
# In <project_dir>
cd src/
nohup bash STAR_pipeline.sh --pair \
-d ~/RNApipeTest \
-i ~/RNApipeTest/data \
-o ~/RNApipeTest/results \
--ref /Reference/mouse/index/star_index/GRCm38.p5_primary_assembly/ \
--gtf /Reference/mouse/annotation/gencode.vM23.annotation.gtf \
-t 16 >nohup1.out 2>&1 &
```

Messages from the program will be direct to file `nohup1.out`, you can use `cat nohup1.out` or `tail nohup1.out` to check. 

 

# Output

1. The trimmed `fastqs` will be output in `<data_dir>/clean` 

```shell
data/clean/
|-- *_trimmed.fq.gz
...
```

2. The QC results of trimmed reads are in ` <output_dir>/QC/`

```shell
results/QC/
|-- *_R1_clean_fastqc.html
|-- *_R1_clean_fastqc.zip

```

Also, `multiqc` reports of `FastQC` will be generated at `results/QC/`. 



3.  The mapping results in `<output_dir>/star/`

```shell
results/star/
|-- *_align
|   |-- Aligned.out.sam
|   |-- Log.final.out
|   |-- Log.out
|   |-- Log.progress.out
|   |-- *.sorted.bam
|   |-- *.sorted.bam.bai
|   -- SJ.out.tab
...

```

Also, `multiqc` reports of `Log.final.out` from `STAR` alignment will be generated at`results/star/`. 



4. The gene expression quantification results in `<output_dir>/featurecounts/`

```shell
  results/featurecounts/
|-- *_counts
|   |-- *_counts.txt
|   `-- *_counts.txt.summary
```

Counts files in `results/featurecounts/*_counts/*_counts.txt` can be used for further analysis. 



# Citations

If you use `STAR_pipeline.sh` for your analysis, please cite the publication: [ONE-seq: epitranscriptome and gene-specific profiling of NAD-capped RNA]()

And the tools used by the pipeline: 

> 1. Krueger, F. (2015) Trim galore. *A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files*, **516**, 517.
> 2. Krueger, F. (2015) Trim galore. *A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files*, **516**, 517.
> 3. Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R. (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, **29**, 15-21.
> 4. Liao, Y., Smyth, G.K. and Shi, W. (2014) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, **30**, 923-930.






- **Author**: Dean Li
- **Date**: 2022-01-26


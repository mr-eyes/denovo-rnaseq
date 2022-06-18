# denovo-rnaseq

## Description
- This is a simple Snakemake workflow to perform de novo transcriptome assembly, DGE, and annotation.
- The workfow currently annotates the highly differential expressed genes with PVal of 0.005 and log_fold_change of 1.

## Create conda enviornment

```
conda env create -f environment.yml
conda activate rnaseq
```

## Limitations
- This workflow is meant to be easy to run, fast to implement, so it's currently not supporting flexibile configs.
- If you need to change configuration and parameters you will need to edit the Snakefile.
- It's only supporting RNASeq paired-end gzipped samples (not interleaved).

## Main software used
- Error trimming: fastp
- De novo transcriptome assembly: rnaSpades
- Quantification: Salmon
- Differential gene expression: DESEQ2
- Annotation: Trinotate

## How to use?

1. Create a working directory
2. Change the `ROOT_DIR` in the `Snakefile` to match your working directory.
3. Create a directory with the name `samples` inside the working directory.
4. Put your samples in the `samples` directory with the naming convection:
    - R1: <sample_name>_1.fastq.gz
    - R2: <sample_name>_2.fastq.gz
5. Copy paste the tab-delimited file [samples.tsv](workflow/samples.tsv) in your workflow directory.
6. Modify the `samples.tsv` to match your samples. Columns as following (sample_type, sample_name, R1_path, R2_path).

## Workflow rulegraph

![](rulegraph.png?raw=true)
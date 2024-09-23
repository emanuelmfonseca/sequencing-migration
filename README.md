# Genomic Sequencing Migration Pipeline

This repository contains a bioinformatics pipeline designed to process simulated genomic sequencing data. The pipeline automates the processing of FASTQ files and manages the data using a cloud infrastructure. It is optimized for scalability and efficiency, reducing the need for manual intervention.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Pipeline Overview](#pipeline-overview)
- [Input Data](#input-data)
- [Usage](#usage)
- [Outputs](#outputs)
- [AWS Architecture](#aws-architecture)
- [Cost Estimate](#cost-estimate)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This pipeline automates the migration and processing of simulated genomic sequencing data, providing robust and scalable solutions for handling large datasets. It is built using Snakemake and designed to run on AWS infrastructure, leveraging cloud services to optimize performance and reduce processing time.

## Installation

### Dependencies

To run this project, you will need the following dependencies installed:

- **Conda**: An environment manager to handle package installation and dependency management.
- **AWS CLI**: Command-line tool for interacting with Amazon Web Services.

### Installation Instructions

You can install all the required dependencies by following these steps:

1. Ensure you have Conda installed. If not, you can download it from [here](https://docs.conda.io/en/latest/miniconda.html).

2. Clone the repository and navigate to the project directory:

    ```bash
    git clone https://github.com/emanuelmfonseca/sequencing-migration.git
    cd sequencing-migration
    ```

3. Create and activate the Conda environment using the provided `environment.yml` file:

    ```bash
    conda env create -f envs/environment.yml
    conda activate sequencing-migration
    ```

This will automatically install **Snakemake** and any other necessary packages specified in the environment configuration. **Snakemake** is a powerful workflow management system designed to create reproducible, scalable, and automated data analysis pipelines.

4. Install the AWS CLI if it's not included in your Conda environment:

    ```bash
    conda install -c conda-forge awscli
    ```

### Setup AWS Credentials
Configure your AWS credentials for S3 access:

```bash
aws configure
```

## Pipeline Overview
The pipeline processes simulated sequencing data by performing primary and secondary analysis stages, including demultiplexing, alignment, and variant calling.

1. ### Primary Analysis

The primary analysis focuses on the initial steps of preparing the data. This phase is crucial because it ensures that we generate high-quality datasets, including trimmed reads and assessments for adapter sequences. By establishing these initial steps, we improve the quality of the subsequent analysis.

- **Generate SNP Matrices**
    - The first step involves running a demographic simulation script to create SNP matrices for multiple weekly runs. This simulates genetic variation based on input demographic parameters.

- **Generate Genomes**
    - In this stage, reference and individual genomes are generated from the SNP matrices. This involves using another simulation script to produce FASTA files for each weekly run, creating both reference genomes and individual genomes for the specified number of samples.

- **Concatenate Individual Genomes**
    - Individual genome FASTA files are concatenated into multi-FASTA files for each weekly run. This prepares the data for generating reads in the next step.

- **Simulate Reads**
-    The multi-FASTA genomes are used to generate simulated FASTQ files (R1 and R2) using the InSilicoSeq tool. This step simulates sequencing reads for the genome.

- **Quality Assessment of Raw Reads**
    - FastQC reports are generated for the simulated FASTQ files to assess the quality of the sequencing data.

- **Trimming Reads**
    - Trimmomatic is utilized to trim the raw FASTQ files to remove low-quality bases and adapter sequences, producing cleaned FASTQ files.

- **Quality Assessment of Trimmed Reads**
    - FastQC reports are generated again for the trimmed reads to evaluate the quality after trimming.

- **Demultiplexing**
    - The trimmed FASTQ files are demultiplexed into sample-specific files based on predefined patterns, allowing for separate analysis of each sample's reads.


2. ### Secondary Analysis

The secondary analysis begins with indexing the reference genome, which organizes it for efficient access. This allows individual sample reads to be mapped to the reference genome, aligning their sequences accurately. Once the reads are aligned, genetic variants can be identified by comparing the sample data to the reference, enabling detailed analysis of genetic differences within the dataset.

- **Index Reference Genome**
    - The reference genome is indexed using BWA to prepare for alignment. This step includes generating necessary auxiliary files (like the FASTA index and sequence dictionary).

- **Align Reads**
    - The demultiplexed FASTQ files are aligned to the reference genome using BWA. The results are stored in SAM and BAM formats, with sorting and duplicate marking performed for each sample.

- **Variant Calling**
    - GATK's HaplotypeCaller is used for variant calling on the aligned BAM files, producing GVCF files for each sample.

- **Combine GVCFs and Genotype**
    - GVCF files from all samples are combined, and GATK's GenotypeGVCFs is run for final variant calling, resulting in a comprehensive VCF file for the weekly run.

## Usage

### Input Data
The input data for this pipeline is automatically generated from simulations of genomic sequencing processes, eliminating the need for external input. However, modifying the configuration file is necessary, particularly to update the working directory.

**Edit Configuration File** 
   Modify the `config.yaml` file to define the necessary paths and key parameters:

   #### Main Parameters
      - working_directory: Specify the full path to the working directory for your project
      - weekly_runs: Number of independent genomic runs to perform
      - nsamples: Number of individuals (samples) processed in each run
      - ncpus: Number of available CPU cores for parallel processing

**Run the Pipeline**:
   To execute the pipeline, run the following command:

```bash
snakemake --cores <number_of_cores>
```

   **Example**:
      If running with 4 cores, use:

```bash
snakemake --cores 4
```

### Output

- The pipeline generates various outputs at each stage, including:

  - **SNP matrix**: A matrix of single nucleotide polymorphisms for each sample, used for population genetic analyses.
    - Example: `data/demographic-simulations/snp-matrix_1.txt`.

  - **Reference genome**: The reference genome used for alignment and variant calling.
    - Example: `data/genome-simulations/genomes_1/fasta/reference_genome_1.fasta`.

  - **Quality reports**: FastQC HTML files for raw and trimmed reads.
    - Example: `data/genome-simulations/genomes_1/qc/genome_1_R1_fastqc.html`.
  
  - **Trimmed FASTQ files**: FASTQ files after trimming.
    - Example: `data/genome-simulations/genomes_1/fastq_trimmed/genome_1_R1_trimmed.fastq`.

  - **Demultiplexed FASTQ files**: FASTQ files split by sample.
    - Example: `data/genome-simulations/genomes_1/demultiplexed/sample_1/genome_1_Ind1_R1_demultiplexed.fastq`.

  - **Aligned SAM and BAM files**: Aligned sequence data and associated metrics.
    - Example: `data/genome-simulations/genomes_1/demultiplexed/sample_1/dedup_sorted_aligned_genome_1.bam`.

  - **Index files**: Index files for efficient BAM file access.
    - Example: `data/genome-simulations/genomes_1/demultiplexed/sample_1/dedup_sorted_aligned_genome_1.bam.bai`.

  - **Variant call files (VCFs)**: Individual and combined VCFs.
    - Example: `data/fully-processed-data/genomes_fully_processed_run_1.vcf`.

## AWS Architecture
The pipeline can be deployed on AWS, utilizing:
- Amazon S3 for data storage
- Amazon EC2 for compute resources
- AWS Batch for job orchestration
- AWS Lambda for automation of processes

An architecture diagram and detailed explanation can be found in the `architecture_diagram.png` file.

## Cost Estimate
A cost estimate for running the pipeline on AWS has been calculated using the AWS Pricing Calculator. This includes:
- EC2 instances for processing
- S3 for data storage
- AWS Batch for workflow orchestration

Details can be found in the `cost_estimate.txt` file.

## Contributing
Contributions are welcome. Please open an issue or submit a pull request for any improvements.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

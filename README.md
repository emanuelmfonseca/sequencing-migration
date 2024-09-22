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
Ensure you have the following dependencies installed:
- Python (>= 3.6)
- Snakemake (>= 7.0)
- AWS CLI
- Conda

You can install the required dependencies by running:

```bash
conda create -n genomics-pipeline python=3.8 snakemake awscli
conda activate genomics-pipeline
```

### Setup AWS Credentials
Configure your AWS credentials for S3 access:

```bash
aws configure
```

### Clone the Repository
Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/yourusername/genomics-pipeline.git
cd genomics-pipeline
```

## Pipeline Overview
The pipeline processes simulated sequencing data by performing primary and secondary analysis stages, including demultiplexing, alignment, and variant calling.

### Stages:
1. **Primary Analysis**: Demultiplexes raw FASTQ files.
2. **Secondary Analysis**: Aligns reads to a reference genome and calls variants.

## Input Data
The input data for this pipeline comes from **simulations** of genomic sequencing processes.

- **Simulated Data**: FASTQ files generated from simulations.
- **Reference Genome**: Reference genome in FASTA format.

### Example input:
- `data/genome-simulations/genomes_1/fastq/genome_1_R1.fastq`
- `data/genome-simulations/genomes_1/fasta/reference_genome_1.fasta`

## Usage
### Step-by-Step Instructions
1. **Edit Configuration File**: 
   Modify the `config.yaml` file to specify the input files and paths:

#### Main Parameters
- working_directory: "/path/to/github_folder"  # Specify the full path to the working directory for your project
- weekly_runs: 5  # Number of independent genomic runs to perform each week
- nsamples: 100  # Number of individuals (samples) processed in each run
- ncpus: 8  # Number of available CPU cores for parallel processing


2. **Run the Pipeline**:
   To execute the pipeline, run the following command:

```bash
snakemake --cores <number_of_cores>
```

3. **Example**:
   If running with 4 cores, use:

```bash
snakemake --cores 4
```

## Outputs
The pipeline generates the following outputs:
- **BAM files**: Aligned sequence data.
- **VCF files**: Variant calls.

### Example output:
- `data/genome-simulations/genomes_1/demultiplexed/sample_1/dedup_sorted_aligned_genome_1.bam`
- `data/fully-processed-data/genomes_fully_processed_run_1.vcf`

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

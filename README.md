# Genomic Sequencing Migration Pipeline

This repository contains a bioinformatics pipeline designed to process simulated genomic sequencing data. The pipeline automates the processing of FASTQ files and manages the data using a cloud infrastructure. It is optimized for scalability and efficiency, reducing the need for manual intervention.

## Table of Contents
- [Introduction](#introduction)
- [Pipeline Overview](#pipeline-overview)
- [AWS Architecture](#aws-architecture)
- [Cost Estimate](#cost-estimate)
- [Installation](#installation)
- [Input Data](#input-data)
- [Usage](#usage)
- [Outputs](#outputs)
- [Development and Testing Environment](#development-and-testing-environment)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This pipeline automates the migration and processing of simulated genomic sequencing data, providing robust and scalable solutions for handling large datasets. It is built using Snakemake and designed to run on AWS infrastructure, leveraging cloud services to optimize performance and reduce processing time.

## Genomic Pipeline Overview
The pipeline processes simulated sequencing data by performing primary and secondary analysis stages, including demultiplexing, alignment, and variant calling.

![genomic-pipeline](https://github.com/user-attachments/assets/09d7889a-4f5b-40c3-94a5-8ac4ed0c7a40)

1. ### Primary Analysis

The primary analysis focuses on the initial steps of preparing the data. This phase is crucial because it ensures that we generate high-quality datasets, including trimmed reads and assessments for adapter sequences. By establishing these initial steps, we improve the quality of the subsequent analysis.

- **Generate SNP Matrices**
    - The first step involves running a demographic simulation script to create SNP matrices for multiple weekly runs. This simulates genetic variation based on input demographic parameters.

- **Generate Genomes**
    - In this stage, reference and individual genomes are generated from the SNP matrices. This involves using another simulation script to produce FASTA files for each weekly run, creating both reference genomes and individual genomes for the specified number of samples.

- **Concatenate Individual Genomes**
    - Individual genome FASTA files are concatenated into multi-FASTA files for each weekly run. This prepares the data for generating reads in the next step.

- **Simulate Reads**
  - The multi-FASTA genomes are used to generate simulated FASTQ files (R1 and R2) using the InSilicoSeq tool. This step simulates sequencing reads for the genome.

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

## AWS Architecture

### **Genomic Data Processing Pipeline Using AWS with IAM and Quilt Integration**:
This genomic pipeline is designed to automate sequencing data processing using AWS services, Quilt for data management, and secure access management through **IAM**. The workflow covers data ingestion, processing, and storage, managed entirely through AWS infrastructure with Quilt providing version control and data lineage.

![aws-architecture-figure](https://github.com/user-attachments/assets/41d5386f-952e-4b98-8b2a-7cfd281f1066)


#### **1. Data Ingestion into S3 and Quilt**:
- Sequencing data is uploaded to an **S3 bucket** and registered in **Quilt** for data management. S3 serves as the central storage for raw sequencing data (e.g., FASTQ files), intermediate results, and final outputs, while Quilt ensures metadata tracking and versioning.
- **IAM roles** control access to S3 and Quilt, allowing only authorized services and users to read and write data.
- Every time new data is uploaded, an event triggers AWS services to start the pipeline.

#### **2. Lambda Trigger**:
- **AWS Lambda** is automatically triggered when new data is uploaded to **S3** and registered in Quilt. Lambda functions, using **IAM roles**, have the permissions required to interact with both **S3** and Quilt, execute workflows, and trigger downstream processes (e.g., EC2 or AWS Batch).
- This automation ensures the pipeline starts processing new data as soon as it becomes available.

#### **3. Running the Snakemake Pipeline on EC2**:
- **Lambda**, using its designated **IAM role**, triggers the launch of an **EC2 instance** to execute the Snakemake pipeline. The **EC2 instance** is assigned an **IAM role** with permissions to access **S3** and **Quilt** for retrieving raw data, running the pipeline, and sending logs and metrics to **CloudWatch**.
- The pipeline, running on EC2, downloads raw data from S3 via Quilt, processes it (e.g., quality control, alignment, variant calling), and tracks metadata updates.

#### **4. Storing Outputs in S3 and Quilt**:
- After processing, output files (e.g., BAM, VCF files) are uploaded back to **S3** and registered in **Quilt** for versioning and data management.
- The **EC2 instance**, using its **IAM role**, ensures that results are securely stored in both **S3** and Quilt for traceability.

#### **5. Batch Processing**:
- For larger datasets, **AWS Step Functions** will be used to orchestrate parallel processing. Step Functions will manage the workflow by submitting jobs to **AWS Batch**, which provisions **EC2** instances to run the Snakemake pipeline. **IAM roles** assigned to the **AWS Batch** instances will ensure secure access to **S3** data, interaction with **Quilt** for metadata management, and the execution of the Snakemake workflows.

#### **6. Workflow Orchestration with Step Functions**:
- **AWS Step Functions** orchestrate the pipeline, ensuring tasks such as data retrieval, processing, and storage are completed sequentially.
- **IAM roles** allow Step Functions to manage permissions between Lambda, EC2, Quilt, and other services, enabling secure and orderly workflow execution.

#### **7. Monitoring and Logging with CloudWatch**:
- Logs from EC2, Lambda, and Batch instances are sent to **CloudWatch** for real-time monitoring.
- **IAM roles** ensure secure transmission of logs and metrics, enabling detailed monitoring and troubleshooting.

#### **8. Secure Access with IAM**:
- **IAM roles** govern access between AWS services and Quilt, ensuring that each service (Lambda, EC2, Batch) has the necessary permissions to interact with other AWS resources securely.
- **IAM roles** are assigned across the pipeline, providing granular control over access to S3, Quilt, Lambda, EC2, Step Functions, and CloudWatch to ensure robust security. 

### Updated **Data Flow Summary with Quilt**:
1. **S3 → Quilt → Lambda**: Data is uploaded to S3 and registered in Quilt. This triggers a Lambda function, with **IAM roles** granting access to both S3 and Quilt.
2. **Lambda → EC2**: Lambda starts an EC2 instance, and **IAM roles** allow the EC2 to fetch data from Quilt (via S3) and run the Snakemake pipeline.
3. **EC2 → Quilt → S3**: Processed results are registered back into Quilt and stored in S3 using the **IAM role** assigned to EC2.
4. **EC2 → CloudWatch**: Logs from EC2 are sent to CloudWatch using **IAM permissions**.
5. **Optional**: **AWS Batch** uses **IAM roles** to manage job scheduling for parallel processing.
6. **Step Functions → Lambda/EC2**: Step Functions orchestrate the workflow using **IAM roles** for managing permissions between all services.

## Cost Estimate

**Here are preliminary cost estimates for the key AWS services involved in the genomic pipeline, as calculated using the AWS Pricing Calculator:**

| Service                    | Details                                                                                                         | Price per Run ($) | Value per Month ($) |
|----------------------------|-----------------------------------------------------------------------------------------------------------------|-------------------|---------------------|
| **Amazon S3 (Storage)**     | - **Storage**: 1.1 TB of data per month <br> - **Price for S3 Standard**: ~$0.023 per GB/month <br> - **Intermediate files**: SAM and BAM typically require 2-3 times more storage than the original raw data <br> - **Monthly storage cost**: $77.72 (3 * 1.1 TB * 1024 GB/TB * $0.023) | 2.43              | 77.72               |
| **Amazon S3 (Retrieval)**   | - **Data Retrieval**: 1.1 TB retrieved per month <br> - **Standard retrieval price**: ~$0.01 per GB <br> - **Retrieval cost**: $11.00 (1.1 TB * 1024 GB/TB * $0.01) | 0.34              | 11.00               |
| **EC2 (Snakemake pipeline)**| - **Instance type**: m5.large (2 vCPUs, 8 GB RAM) <br> - **Estimated hours per sequencing run**: 24 hours per run <br> - **Runs per week**: 8 runs (total 32 runs/month) <br> - **Total hours per month**: 24 hours per run * 32 runs/month = 768 hours/month <br> - **Price for m5.large**: ~$0.096 per hour <br> - **Monthly EC2 cost**: 768 hours * $0.096 = $73.73 | 2.30              | 73.73               |
| **Lambda (Trigger function)**| - **Number of invocations**: 32 per month (1 per sequencing run) <br> - **Execution time**: ~1 second (1000 ms) <br> - **Memory**: 128 MB <br> - **Price**: $0.00001667 per request (100 ms increments) <br> - **Monthly Lambda cost**: $0.05 (32 invocations * $0.00001667) | 0.0016            | 0.05                |
| **AWS Batch**              | - **Assume 10 EC2 instances (m5.large)** running in parallel for batch jobs <br> - **Total runtime for batch jobs**: 24 hours per run <br> - **Total EC2 hours for Batch**: 10 instances * 24 hours = 240 hours per run <br> - **For 32 runs per month**: 240 hours per run * 32 runs = 7,680 hours per month <br> - **Price for m5.large**: ~$0.096 per hour <br> - **Monthly Batch cost**: 7,680 hours * $0.096 = $737.28 | 23.04             | 737.28              |
| **CloudWatch**             | - **Logs volume**: 10 GB/month <br> - **Custom metrics**: 5 metrics/month <br> - **Log storage cost**: ~$0.50 per GB <br> - **Metrics cost**: ~$0.30 per metric <br> - **Monthly CloudWatch cost**: $8.50 (10 GB * $0.50 + 5 metrics * $0.30) | 0.27              | 8.50                |
| **Step Functions**         | - **State transitions**: 32 transitions (1 per sequencing run) <br> - **Cost per transition**: $0.025 per 1,000 transitions <br> - **Monthly Step Functions cost**: $0.01 (32 transitions * $0.025/1000) | 0.00031           | 0.01                |

**Total Cost per Run**: $28.05
**Total Monthly Cost**: $897.29

## Installation

### Dependencies

To run this project, you will need the following dependencies installed:

- **Conda**: An environment manager to handle package installation and dependency management.

### Installation Instructions

**Follow these steps to install Conda and access the pipeline:**

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

## Usage

### Input Data
The input data for this pipeline is automatically generated from simulations of genomic sequencing processes, eliminating the need for external input. However, modifying the configuration file is necessary, particularly to update the working directory.

**Edit Configuration File**: Modify the `config.yaml` file to define the necessary paths and key parameters:

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

### Outputs

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

## Development and Testing Environment
This pipeline was developed on a **MacBook 2020 with an M1 chip**. The tutorial associated with this project was tested on another machine with the same configuration.

## Contributing
Contributions are welcome. Please open an issue or submit a pull request for any improvements. For questions, feel free to email emanuelmfonseca@gmail.com.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

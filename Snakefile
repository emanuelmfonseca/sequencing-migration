configfile: "config.yaml"


rule all:
    input:
        fastqc=expand(
            "data/genome-simulations/genomes_{weekly_run}/qc/genome_{weekly_run}_{read}_fastqc.html",
            weekly_run=range(1, config["weekly_runs"] + 1),
            read=["R1", "R2"]
        ),
        
        fastqc_trimmed=expand(
            "data/genome-simulations/genomes_{weekly_run}/qc_trimmed/genome_{weekly_run}_{read}_trimmed_fastqc.html",
            weekly_run=range(1, config["weekly_runs"] + 1),
            read=["R1", "R2"]),
        
        # Define the paths for the trimmed FASTQ files using expand.
        trimmed_fastq=expand(
            "data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/genome_{weekly_run}_{read}_trimmed.fastq",
            weekly_run=range(1, config["weekly_runs"] + 1),
            read=["R1", "R2"]
        ),
        
        index=expand(
            "data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta.{suffix}",
            weekly_run=range(1, config["weekly_runs"] + 1),
            suffix=['amb', 'ann', 'bwt', 'pac', 'sa', 'fai']
        ),

        index_bam=expand(
            "data/genome-simulations/genomes_{weekly_run}/aligned_genome/sorted_aligned_genome_{weekly_run}.bam.bai",
            weekly_run=range(1, config["weekly_runs"] + 1)
        ),

        dicti=expand(
            "data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.dict",
            weekly_run=range(1, config["weekly_runs"] + 1)
        ),

        vcf=expand(
            "data/genome-simulations/genomes_{weekly_run}/gatk/genome_{weekly_run}.vcf",
            weekly_run=range(1, config["weekly_runs"] + 1)
            )


# Execute the demographic simulation to produce a SNP matrix for multiple runs
rule generate_snps:
    input:
        script=config["working_directory"] + "/scripts/demographic-simulation.py"  # Path to the simulation script
    output:
        expand("data/demographic-simulations/snp-matrix_{rep}.txt", rep=range(1, config["weekly_runs"] + 1))  # Output SNP matrices for each weekly run
    params:
        nref = config["nref"],                # Effective population size in the past
        nu1 = config["nu1"],                  # Current population size relative to the reference population
        t1 = config["t1"],                    # Time relative to the reference population
        nsamples = config["nsamples"],          # Number of individuals to sample
        ploidy = config["ploidy"],            # Ploidy level (2 for diploid)
        mu = config["mu"],                    # Mutation rate
        genome_size = config["genome_size"],  # Length of the genome
        recomb = config["recomb"],            # Recombination rate
        seed = config["seed"]                 # Random seed for reproducibility
    conda:
        # Specify the conda environment for dependencies
        "envs/environment.yml"
    shell:
        "python {input.script} {params.nref} {params.nu1} {params.t1} {params.nsamples} {params.ploidy} {params.mu} {params.genome_size} {params.recomb} {params.seed} {output}"

# Generates reference and individual genomes from SNP matrices
rule generate_genome:
    input:
        # Input SNP matrices from demographic_simulation rule
        snp_matrices=expand("data/demographic-simulations/snp-matrix_{weekly_run}.txt", weekly_run=range(1, config["weekly_runs"] + 1)),  
        # Path to the genome simulation script
        script=config["working_directory"] + "/scripts/genome-simulation.py",
        # Path to the working directory
        working_directory=config["working_directory"]
    output:
        # Output for the reference genome
        reference=expand("data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta", weekly_run=range(1, config["weekly_runs"] + 1)),  
        # Output for individual genomes
        individual=expand("data/genome-simulations/genomes_{weekly_run}/fasta/individual_genome_{ind}.fasta", weekly_run=range(1, config["weekly_runs"] + 1), ind=range(1, config["nsamples"] + 1))  
    conda:
        # Specify the conda environment for dependencies
        "envs/environment.yml" 
    shell:
        # Call the genome simulation script with input and output
        "python {input.script} {input.working_directory} {input.snp_matrices}"


# Rule to concatenate individual genome FASTA files into multi-FASTA files per weekly run.
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Individual genome FASTA files for each run and sample
            individual=expand(f"data/genome-simulations/genomes_{weekly_run}/fasta/individual_genome_{{ind}}.fasta", ind=range(1, config["nsamples"] + 1))
        output:
            # Multi-FASTA output file for each weekly run
            multifasta=f"data/genome-simulations/genomes_{weekly_run}/multifasta/multifasta_genome_{weekly_run}.fasta"
        conda:
            # Specify the conda environment for dependencies
            "envs/environment.yml"  
        shell:
            """
            touch {output}  # Ensure the output file exists
            cat {input.individual} >> {output}  # Concatenate individual genomes
            """


# Rule to generate simulated FASTQ files from multi-FASTA genome inputs for each weekly run.
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Specify the multi-FASTA file that contains genome sequences for the current weekly run.
            multifasta=f"data/genome-simulations/genomes_{weekly_run}/multifasta/multifasta_genome_{weekly_run}.fasta"
        output:
            # Define the paths for the generated FASTQ files (R1 and R2) for each weekly run.
            expand("data/genome-simulations/genomes_{weekly_run}/fastq/genome_{weekly_run}_{read}.fastq", weekly_run=weekly_run, read=["R1", "R2"])
        params:
            # Store the current weekly run number and number of CPUs as a parameter for use in the shell commands.
            weekly_run=weekly_run,
            ncpus=config["ncpus"]
        conda:
            # Specify the Conda environment to ensure that all dependencies for the rule are installed.
            "envs/environment.yml"
        shell:
            """
            # Generate simulated reads using InSilicoSeq (iss) with the HiSeq model.
            iss generate --model hiseq --genomes {input.multifasta} --n_reads 1000 --cpus {params.ncpus} --output data/genome-simulations/genomes_{params.weekly_run}/fastq/genome_{params.weekly_run}
            
            #Remove any unnecessary .txt and .vcf files generated by iss.
            rm data/genome-simulations/genomes_{params.weekly_run}/fastq/*.txt data/genome-simulations/genomes_{params.weekly_run}/fastq/*.vcf
            """


# Rule to generate FastQC reports for quality assessment of sequencing data in each weekly run.
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Specify the paths for the input FASTQ files (R1 and R2) for each weekly run.
            fastq=expand("data/genome-simulations/genomes_{weekly_run}/fastq/genome_{weekly_run}_{read}.fastq", weekly_run=weekly_run, read=["R1", "R2"])
        output:
            # Define the paths for the FastQC output HTML files (one for each read).
            expand("data/genome-simulations/genomes_{weekly_run}/qc/genome_{weekly_run}_{read}_fastqc.html", weekly_run=weekly_run, read=["R1", "R2"])
        params:
            # Store the current weekly run number for use in shell commands.
            weekly_run=weekly_run,
        conda:
            # Specify the Conda environment.
            "envs/environment.yml"
        shell:
            # Run FastQC on the input files, save output in the qc folder for the current run.
            """
            fastqc {input} -o data/genome-simulations/genomes_{params.weekly_run}/qc
            """

# Rule to run Trimmomatic for trimming FASTQ files in each weekly run.
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Specify the input FASTQ files (R1 and R2) for trimming.
            fastq=[
                f"data/genome-simulations/genomes_{weekly_run}/fastq/genome_{weekly_run}_R1.fastq",
                f"data/genome-simulations/genomes_{weekly_run}/fastq/genome_{weekly_run}_R2.fastq"
            ]
        output:
            # Define the paths for the trimmed FASTQ files.
            trimmed_fastq=[
                f"data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/genome_{weekly_run}_R1_trimmed.fastq",
                f"data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/genome_{weekly_run}_R2_trimmed.fastq"
            ],
            # Define the path for the log file with the same wildcard as weekly_run.
            log=f"data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/trimmomatic_{weekly_run}.log"
        params:
            # Parameters for Trimmomatic (modify based on your needs).
            trimmomatic_params="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            weekly_run=weekly_run
        conda:
            # Specify the Conda environment for Trimmomatic.
            "envs/environment.yml"
        shell:
            # Run Trimmomatic to trim the input FASTQ files.
            """
            trimmomatic PE -phred33 {input.fastq[0]} {input.fastq[1]} \
            {output.trimmed_fastq[0]} {output.trimmed_fastq[0]}.unpaired \
            {output.trimmed_fastq[1]} {output.trimmed_fastq[1]}.unpaired \
            {params.trimmomatic_params} \
            > {output.log} 2>&1
            """

# Rule to generate FastQC reports from trimmed FASTQ files for quality reassessment
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Specify the paths for the input FASTQ files (R1 and R2) for each weekly run.
            fastq_trimmed=expand("data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/genome_{weekly_run}_{read}_trimmed.fastq", weekly_run=weekly_run, read=["R1", "R2"])
        output:
            # Define the paths for the FastQC output HTML files (one for each read).
            expand("data/genome-simulations/genomes_{weekly_run}/qc_trimmed/genome_{weekly_run}_{read}_trimmed_fastqc.html", weekly_run=weekly_run, read=["R1", "R2"])
        params:
            # Store the current weekly run number for use in shell commands.
            weekly_run=weekly_run,
        conda:
            # Specify the Conda environment.
            "envs/environment.yml"
        shell:
            # Run FastQC on the input files, save output in the qc folder for the current run.
            """
            fastqc {input} -o data/genome-simulations/genomes_{params.weekly_run}/qc_trimmed
            """


# Align trimmed FASTQ files using BWA for each weekly run
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Reference genome and paired-end FASTQ files (R1, R2)
            reference=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta",
            fastq_trimmed=expand("data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/genome_{weekly_run}_{read}_trimmed.fastq", weekly_run=weekly_run, read=["R1", "R2"])
        output:
            # Indexed reference files and sorted BAM
            index=expand("data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta.{suffix}", weekly_run=weekly_run, suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
            fai=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta.fai",
            dicti=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.dict",
            sam=f"data/genome-simulations/genomes_{weekly_run}/aligned_genome/aligned_genome_{weekly_run}.sam",
            bam=f"data/genome-simulations/genomes_{weekly_run}/aligned_genome/aligned_genome_{weekly_run}.bam",
            sorted_bam=f"data/genome-simulations/genomes_{weekly_run}/aligned_genome/sorted_aligned_genome_{weekly_run}.bam",
            index_bam=f"data/genome-simulations/genomes_{weekly_run}/aligned_genome/sorted_aligned_genome_{weekly_run}.bam.bai"
        params:
            weekly_run=weekly_run,
            
            # Define the read group parameters
            RGID=f"tsk_{weekly_run}",    # Read group identifier
            RGLB=f"tsk_{weekly_run}",    # Library identifier
            RGPL="ILLUMINA",             # Platform
            RGPM="HISEQ",                # Platform model
            RGSM=f"tsk_{weekly_run}"     # Sample name

        conda:
            # Use BWA in Conda environment
            "envs/environment.yml"
        threads: config["ncpus"]  # Set number of threads for BWA
        shell:
            """
            # Index the reference genome
            bwa index {input.reference}
            
            # Generate the FASTA index (.fai) file using samtools faidx
            samtools faidx {input.reference}
            
            # Generate the sequence dictionary (.dict) file using Picard
            picard CreateSequenceDictionary R={input.reference} O={output.dicti}
            
            # Align the trimmed FASTQ reads to the reference genome and output SAM file
            bwa mem -M -R "@RG\\tID:{params.RGID}\\tLB:{params.RGLB}\\tPL:{params.RGPL}\\tPM:{params.RGPM}\\tSM:{params.RGSM}" {input.reference} {input.fastq_trimmed} > {output.sam}
            
            # Convert SAM to BAM
            samtools view -Sb {output.sam} > {output.bam}
            
            # Sort the BAM file
            samtools sort {output.bam} -o {output.sorted_bam}

            # Index the sorted BAM file
            samtools index {output.sorted_bam}
            """


# Rule to run GATK HaplotypeCaller on sorted BAM files for variant calling
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Reference genome and sorted BAM file from the previous step
            reference=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta",
            bam=f"data/genome-simulations/genomes_{weekly_run}/aligned_genome/sorted_aligned_genome_{weekly_run}.bam",
        output:
            # GATK output: VCF file containing variants
            vcf=f"data/genome-simulations/genomes_{weekly_run}/gatk/genome_{weekly_run}.vcf"
        conda:
            # Use GATK in Conda environment
            "envs/environment.yml"
        threads: config["ncpus"]  # Set number of threads for GATK
        shell:
            """
            # Run GATK HaplotypeCaller for variant calling
            gatk HaplotypeCaller --reference {input.reference} --input {input.bam} --output {output.vcf}
            """



# Run GATK HaplotypeCaller for variant calling
#gatk HaplotypeCaller --reference data/genome-simulations/genomes_1/fasta/reference_genome_1.fasta --input data/genome-simulations/genomes_1/aligned_genome/sorted_aligned_genome_1.bam --output data/genome-simulations/genomes_1/gatk/genome_1.vcf

#picard CreateSequenceDictionary R=/Users/emanuelmfonseca/project/sequencing-migration/data/genome-simulations/genomes_1/fasta/reference_genome_1.fasta O=data/genome-simulations/genomes_1/fasta/reference_genome_1.fasta.dict

#bwa mem mem -M -R @RG\tID:tsk_1 \tLB:tsk_1 \tPL:ILLUMINA\tPM:HISEQ\tSM:tsk_1 data/genome-simulations/genomes_2/fasta/reference_genome_2.fasta data/genome-simulations/genomes_2/fastq_trimmed/genome_2_R1_trimmed.fastq data/genome-simulations/genomes_2/fastq_trimmed/genome_2_R2_trimmed.fastq > data/genome-simulations/genomes_2/aligned_genome/aligned_genome_2.sam
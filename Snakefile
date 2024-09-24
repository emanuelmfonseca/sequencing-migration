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

        demu=expand(
            "data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/genome_{weekly_run}_Ind{sample}_{read}_demultiplexed.fastq",
            weekly_run=range(1, config["weekly_runs"] + 1),
            sample=range(1, config['nsamples']+1),
            read=["R1", "R2"]),
        
        index=expand(
            "data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta.{suffix}",
            weekly_run=range(1, config["weekly_runs"] + 1),
            suffix=['amb', 'ann', 'bwt', 'pac', 'sa', 'fai']
        ),
                
        dicti=expand(
            "data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.dict",
            weekly_run=range(1, config["weekly_runs"] + 1)
        ),

        index_bam=expand(
            "data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/dedup_sorted_aligned_genome_{sample}.bam.bai",
            weekly_run=range(1, config["weekly_runs"] + 1),
            sample=range(1, config["nsamples"] + 1)
        ),

        g_vcf=expand(
            "data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/sample_{weekly_run}_{sample}.g.vcf",
            weekly_run=range(1, config["weekly_runs"] + 1),
            sample=range(1, config["nsamples"] + 1)
        ),

        final_vcf=expand(
            "data/fully-processed-data/genomes_fully_processed_run_{weekly_run}.vcf",
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
        "environment.yml"
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
        "environment.yml" 
    shell:
        # Call the genome simulation script with input and output
        "python {input.script} {input.working_directory} {input.snp_matrices}"


# Rule to concatenate individual genome FASTA files into multi-FASTA files per weekly run.
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Individual genome FASTA files for each run and sample
            individual=expand(
                "data/genome-simulations/genomes_{weekly_run}/fasta/individual_genome_{ind}.fasta",
                weekly_run=weekly_run,
                ind=range(1, config["nsamples"] + 1))
        output:
            # Multi-FASTA output file for each weekly run
            multifasta=f"data/genome-simulations/genomes_{weekly_run}/multifasta/multifasta_genome_{weekly_run}.fasta"
        params:
            nsamples=config["nsamples"]
        conda:
            # Specify the conda environment for dependencies
            "environment.yml"  
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
            ncpus=config["ncpus"],
            nreads=int(float(config["nreads"])),
        conda:
            # Specify the Conda environment to ensure that all dependencies for the rule are installed.
            "environment.yml"
        shell:
            """
            # Generate simulated reads using InSilicoSeq (iss) with the HiSeq model.
            iss generate --model hiseq --genomes {input.multifasta} --n_reads {params.nreads} --cpus {params.ncpus} --output data/genome-simulations/genomes_{params.weekly_run}/fastq/genome_{params.weekly_run}
            
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
            "environment.yml"
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
            "environment.yml"
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
            "environment.yml"
        shell:
            # Run FastQC on the input files, save output in the qc folder for the current run.
            """
            fastqc {input} -o data/genome-simulations/genomes_{params.weekly_run}/qc_trimmed
            """


# This rule demultiplexes trimmed FASTQ files by sample for each weekly run.
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Specify the paths for the input FASTQ files (R1 and R2) for each weekly run.
            fastq_trimmed=expand(
                "data/genome-simulations/genomes_{weekly_run}/fastq_trimmed/genome_{weekly_run}_{read}_trimmed.fastq",
                weekly_run=weekly_run, 
                read=["R1", "R2"]
            )
        output:
            # Define the paths for the demultiplexed FASTQ files (one for each pattern and read).
            demultiplexed=expand(
                "data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/genome_{weekly_run}_Ind{sample}_{read}_demultiplexed.fastq",
                weekly_run=weekly_run,
                sample=range(1, config['nsamples']+1),
                read=["R1", "R2"]
            )
        params:
            # Store the current weekly run number.
            weekly_run=weekly_run,
            # Create patterns based on the number of samples.
            patterns=["Ind" + str(i) for i in range(1, config['nsamples'] + 1)]  # List of patterns, not a string
        conda:
            # Specify the Conda environment.
            "environment.yml"
        shell:
            # Loop over patterns and demultiplex each one using seqkit.
            """
            for pattern in {params.patterns}; do
                grep -A 3 "$pattern" {input.fastq_trimmed[0]} > data/genome-simulations/genomes_{params.weekly_run}/demultiplexed/sample_${{pattern:3}}/genome_{params.weekly_run}_${{pattern}}_R1_demultiplexed.fastq;
                
                grep -A 3 "$pattern" {input.fastq_trimmed[1]} > data/genome-simulations/genomes_{params.weekly_run}/demultiplexed/sample_${{pattern:3}}/genome_{params.weekly_run}_${{pattern}}_R2_demultiplexed.fastq;
            done
            """

# Index reference genome and create auxiliary files for each weekly run
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Reference genome
            reference=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta",
        output:
            # Indexed reference files
            index=expand("data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta.{suffix}", weekly_run=weekly_run, suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
            fai=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta.fai",
            dicti=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.dict",
        conda:
            # Specify the Conda environment.
            "environment.yml"
        shell:
            """
            # Index the reference genome
            bwa index {input.reference}
            
            # Generate the FASTA index (.fai) file using samtools faidx
            samtools faidx {input.reference}
            
            # Generate the sequence dictionary (.dict) file using Picard
            picard CreateSequenceDictionary R={input.reference} O={output.dicti}
            """


# Align trimmed FASTQ files using BWA for each weekly run
for weekly_run in range(1, config["weekly_runs"] + 1):
    for sample in range(1, config["nsamples"] + 1):
        rule:
            input:
                # Reference genome file and paired-end FASTQ files (R1 and R2)
                reference=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta",
                
                # Expand file paths for demultiplexed FASTQ files (R1, R2) for each sample
                fastq_demu=expand("data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/genome_{weekly_run}_Ind{sample}_{read}_demultiplexed.fastq", weekly_run=weekly_run, sample=sample, read=["R1", "R2"])
            output:
                # Output SAM, BAM, sorted BAM, and BAM index files
                sam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/aligned_genome_{sample}.sam",
                bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/aligned_genome_{sample}.bam",
                sorted_bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/sorted_aligned_genome_{sample}.bam",
                metrics=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/dedup_metrics.txt",
                dedup_sorted_bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/dedup_sorted_aligned_genome_{sample}.bam",
                index_bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/dedup_sorted_aligned_genome_{sample}.bam.bai"
            params:
                # Parameters for each sample's read group in BWA
                RGID=f"tsk_{sample}",  # Read group identifier
                RGLB=f"tsk_{sample}",  # Library identifier
                RGPL="ILLUMINA",       # Sequencing platform
                RGPM="HISEQ",          # Platform model
                RGSM=f"tsk_{sample}"   # Sample name
            conda:
                # Specify the Conda environment.
                "environment.yml"
            threads: config["ncpus"]  # Number of threads for BWA
            shell:
                """
                # Align FASTQ reads to the reference genome, output SAM file
                bwa mem -M -R "@RG\\tID:{params.RGID}\\tLB:{params.RGLB}\\tPL:{params.RGPL}\\tPM:{params.RGPM}\\tSM:{params.RGSM}" {input.reference} {input.fastq_demu} > {output.sam}
                
                # Convert SAM to BAM
                samtools view -Sb {output.sam} > {output.bam}

                # Sort the BAM file
                samtools sort {output.bam} -o {output.sorted_bam}
                
                picard MarkDuplicates I={output.sorted_bam} O={output.dedup_sorted_bam} M={output.metrics} REMOVE_DUPLICATES=true
                
                # Index the sorted BAM file
                samtools index {output.dedup_sorted_bam}
                """

# Rule to run GATK HaplotypeCaller for variant calling on each sample across multiple weekly runs.
for weekly_run in range(1, config["weekly_runs"] + 1):
    for sample in range(1, config["nsamples"] + 1):
        rule:
            input:
                # Reference genome and sorted BAM file for the current sample.
                reference=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta",
                bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/dedup_sorted_aligned_genome_{sample}.bam",
                index_bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/dedup_sorted_aligned_genome_{sample}.bam.bai"
            output:
                # Output VCF file with variants for the current sample.
                vcf=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_{sample}/sample_{weekly_run}_{sample}.g.vcf"
            conda:
                # Specify the Conda environment.
                "environment.yml"
            threads: config["ncpus"]  # Number of threads for GATK.
            shell:
                """
                gatk HaplotypeCaller --reference {input.reference} --input {input.bam} --output {output.vcf} -ERC GVCF
                """


# Rule to combine GVCFs and run GATK GenotypeGVCFs for variant calling on combined genomes
for weekly_run in range(1, config["weekly_runs"] + 1):
    rule:
        input:
            # Reference genome file for the current run
            reference=f"data/genome-simulations/genomes_{weekly_run}/fasta/reference_genome_{weekly_run}.fasta",
            # Index BAM file and vcf for validation
            index_bam=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_1/dedup_sorted_aligned_genome_1.bam.bai",
            vcf=f"data/genome-simulations/genomes_{weekly_run}/demultiplexed/sample_1/sample_{weekly_run}_1.g.vcf"
        output:
            # Combined VCF file from multiple GVCFs
            comb_vcf=f"data/genome-simulations/genomes_{weekly_run}/gatk/combined_genomes_run_{weekly_run}.g.vcf",
            # Final processed VCF file
            final_vcf=f"data/fully-processed-data/genomes_fully_processed_run_{weekly_run}.vcf"
        params:
            # Weekly run iteration for file paths
            weekly_run=weekly_run
        conda:
            # Specify the Conda environment.
            "environment.yml"
        threads: config["ncpus"]  # Number of threads for GATK
        shell:
            """
            # Collect all GVCF files for the current run
            bam_files=$(find data/genome-simulations/genomes_{params.weekly_run}/demultiplexed -type f -name "*g.vcf" | sed 's/^/--variant /')
            
            # Combine GVCF files into a single VCF
            gatk CombineGVCFs --reference {input.reference} $bam_files --output {output.comb_vcf}
            
            # Genotype the combined VCF
            gatk GenotypeGVCFs --reference {input.reference} --variant {output.comb_vcf} --output {output.final_vcf}
            """

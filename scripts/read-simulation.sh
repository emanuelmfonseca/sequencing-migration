#!/bin/bash

# Arguments from Snakemake
working_directory_path=$1  # Pass the working directory from the input
weekly_run=$2              # Pass weekly_run from wildcards
ind=$3                     # Pass ind from wildcards

# Define the path where the individual genomes are located
path1="${working_directory_path}/data/genome-simulations/genomes_${weekly_run}/fasta"

# Find all fasta files corresponding to individual genomes
fasta_files=($(find "$path1" -type f -name "individual_genome_${ind}.fasta"))

# Iterate over each fasta file and run the read simulation
for fasta_file in ${fasta_files[@]}
do
  iss generate --model hiseq --genomes "$fasta_file" --n_reads 10000 --cpus 4 --output "${working_directory_path}/data/genome-simulations/genomes_${weekly_run}/fastqc/individual_genome_${ind}"
done

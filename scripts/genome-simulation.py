# Import necessary libraries
import os
import sys
import numpy as np
import pandas as pd

# List of output files to be generated

def main(working_directory, input_files):
    
    # Define base pair probabilities for a human-like genome
    base_probabilities = {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2}
    
    for i, input_file in enumerate(input_files):
        # Path to the simulated SNP data
        simulation_file = os.path.join(input_file)
        
        # Load the SNP simulation data
        simulation_data = pd.read_csv(simulation_file, sep='\s+', header=None)
        
        # Extract SNP positions, rounding them to the nearest integer and converting to numpy array
        snp_positions = np.round(simulation_data.iloc[:, 0]).astype(int).to_numpy()
        
        # Extract alleles for each individual
        individual_alleles = simulation_data.iloc[:, 1:].transpose()
        
        # Create indexes for diploid individuals (two alleles per individual)
        allele1_index = np.arange(0, individual_alleles.shape[0], 2)
        allele2_index = np.arange(1, individual_alleles.shape[0], 2)
        
        # Number of diploid individuals
        num_individuals = individual_alleles.shape[0] // 2
        
        # Generate reference and alternative alleles for each SNP position
        np.random.seed(42)  # Ensure reproducibility
        ref_alleles = np.random.choice(list(base_probabilities), len(snp_positions), replace=True, p=list(base_probabilities.values()))

        # Generate random mutations by selecting an alternate allele that differs from the reference allele.
        # This process assumes a Jukes-Cantor mutation model, which posits equal probabilities for all mutations
        alt_alleles = []
        for ii in range(len(ref_alleles)):
            ref = ref_alleles[ii]
            possible_alt_bases = [b for b in base_probabilities if b != ref]  # Exclude the reference base
            np.random.seed(ii)
            alt_alleles.append(np.random.choice(possible_alt_bases, 1)[0])
        
        # Generate a random reference genome sequence
        np.random.seed(42)
        genome_length = max(snp_positions) + 500  # Add buffer length for genome
        reference_genome = np.random.choice(list(base_probabilities), genome_length, replace=True, p=list(base_probabilities.values()))
        reference_genome[snp_positions] = ref_alleles  # Insert SNPs into the reference genome
        
        # Split reference genome into lines of 70 bases for FASTA format
        reference_genome_str = ''.join(reference_genome)
        reference_genome_lines = [reference_genome_str[i:i + 70] for i in range(0, len(reference_genome_str), 70)]
        
        # Ensure output directory exists
        output_dir = os.path.join(working_directory, f'data/genome-simulations/genomes_{i+1}/fasta')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # Write the reference genome to a FASTA file
        reference_genome_path = os.path.join(output_dir, f'reference_genome_{i+1}.fasta')
        with open(reference_genome_path, 'w') as ref_file:
            _ = ref_file.write(">Reference\n")
            for line in reference_genome_lines:
                _ = ref_file.write(line + '\n')
        
        # Generate individual genomes from msprime simulations
        for individual_idx in range(num_individuals):
            allele1 = individual_alleles.iloc[allele1_index[individual_idx]]
            allele2 = individual_alleles.iloc[allele2_index[individual_idx]]
            genotype = allele1 + allele2
            
            # Generate random genomes for each individual
            np.random.seed(42)
            genome1 = np.random.choice(list(base_probabilities), genome_length, replace=True, p=list(base_probabilities.values()))
            genome2 = genome1.copy()  # Copy for the second allele
            
            # Insert SNPs into the individual genomes based on genotype
            genome1_snps = []
            genome2_snps = []
            for snp_idx in range(len(genotype)):
                if genotype[snp_idx] == 0:
                    # Homozygous reference
                    genome1_snps.append(ref_alleles[snp_idx])
                    genome2_snps.append(ref_alleles[snp_idx])
                elif genotype[snp_idx] == 1:
                    # Heterozygous
                    selected_bases = np.random.choice([ref_alleles[snp_idx], alt_alleles[snp_idx]], size=2, replace=False)
                    genome1_snps.append(selected_bases[0])
                    genome2_snps.append(selected_bases[1])
                else:
                    # Homozygous alternate
                    genome1_snps.append(alt_alleles[snp_idx])
                    genome2_snps.append(alt_alleles[snp_idx])
            
            # Insert SNPs into genome sequence
            genome1[snp_positions] = genome1_snps
            genome2[snp_positions] = genome2_snps
            
            # Prepare FASTA format for individual genomes
            genome1_str = ''.join(genome1)
            genome1_lines = [genome1_str[i:i + 70] for i in range(0, len(genome1_str), 70)]
            
            genome2_str = ''.join(genome2)
            genome2_lines = [genome2_str[i:i + 70] for i in range(0, len(genome2_str), 70)]
            
            # Save individual genomes to FASTA file
            individual_genome_path = os.path.join(output_dir, f'individual_genome_{individual_idx+1}.fasta')
            with open(individual_genome_path, 'w') as ind_file:
                _ = ind_file.write(f">Ind{individual_idx + 1}_1\n")
                for line in genome1_lines:
                    _ = ind_file.write(line + '\n')
                
                _ = ind_file.write(f">Ind{individual_idx + 1}_2\n")
                for line in genome2_lines:
                    _ = ind_file.write(line + '\n')


if __name__ == "__main__":
    # Pass output file paths from command line arguments
    working_directory = sys.argv[1]
    input_files = sys.argv[2:]
    main(working_directory, input_files)

























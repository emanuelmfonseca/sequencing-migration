# Import necessary libraries
import msprime
import math
import os
import sys

def main(nref, nu1, t1, nsample, ploidy, mu, genome_size, recomb, seed, output_files):    
    # Calculate time in generations since the population size change
    time_pop_change = 2 * nref * t1
    ncurr = nref * nu1
    
    # Calculate growth rate
    r_curr = math.log(ncurr / nref) / time_pop_change
    
    # Create a demography model
    demography = msprime.Demography()
    demography.add_population(name="pop1", initial_size=ncurr, growth_rate=r_curr)
    demography.add_population_parameters_change(time=time_pop_change, initial_size=nref, growth_rate=0)
    
    for output_file in output_files:

        # Create output directory if it does not exist
        output_dir = os.path.dirname(output_file)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # Simulate ancestry
        ts = msprime.sim_ancestry(samples={"pop1": nsamples}, ploidy=ploidy, 
                                   demography=demography, sequence_length=genome_size, 
                                   recombination_rate=recomb, random_seed=seed)
        
        # Simulate mutations
        mts = msprime.sim_mutations(ts, rate=mu, discrete_genome=False)
        
        # Write the SNP matrix to the output file
        with open(output_file, 'w') as output:
            for var in mts.variants():
                alleles = var.genotypes
                position = round(var.site.position)
                alleles = [str(i) for i in alleles]
                v = [str(position), '\t'.join(alleles)]
                output.write('\t'.join(v) + '\n')


if __name__ == "__main__":
    nref = float(sys.argv[1])
    nu1 = float(sys.argv[2])
    t1 = float(sys.argv[3])
    nsamples = float(sys.argv[4])
    ploidy = float(sys.argv[5])
    mu = float(sys.argv[6])
    genome_size = float(sys.argv[7])
    recomb = float(sys.argv[8])
    seed = float(sys.argv[9]) 
    output_files = sys.argv[10:]  # Get all output file paths from command-line arguments
    main(nref, nu1, t1, nsamples, ploidy, mu, genome_size, recomb, seed, output_files) # Call the main function with the specified output files


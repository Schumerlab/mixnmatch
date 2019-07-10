#!/bin/sh

#SBATCH --ntasks=1 
#SBATCH --mem=32000
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --job-name=sim-hybrid-genomes
#SBATCH --mail-user=youremail@stanford.edu

module load R
module load perl
module load biology py-biopython
module load boost
module load gsl
export PATH="/home/groups/schumer/shared_bin/lib:$PATH"
export PATH="/home/groups/schumer/shared_bin:$PATH"
export LD_LIBRARY_PATH="/home/groups/schumer/shared_bin/lib:$LD_LIBRARY_PATH"

perl /scratch/groups/schumer/molly/Simulate_Ancestry_HMM/Simulate_hybrid_genomes/simulate_admixed_genomes_v6.pl hybrid_simulation_configuration.cfg

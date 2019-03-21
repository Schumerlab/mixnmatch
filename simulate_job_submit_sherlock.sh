#!/bin/sh

#SBATCH --ntasks=1 
#SBATCH --mem=32000
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --job-name=sim-hybrid-genomes
#SBATCH --mail-user=youremail@stanford.edu

module load R
module load perl

perl simulate_admixed_genomes.pl hybrid_simulation_configuration.cfg

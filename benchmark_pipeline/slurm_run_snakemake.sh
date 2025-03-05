#!/bin/bash
#SBATCH --job-name=benchmark_pipeline
#SBATCH --partition=cpu
#SBATCH --time=01-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G  
#SBATCH --ntasks=1
#SBATCH --output=log/benchmark_pipeline_%j.log

source activate benchmark_llms
snakemake --snakefile Snakefile --profile slurm --rerun-incomplete

#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --partition exacloud
#SBATCH --time=10:59:59

ENV_NAME="HOP"

if conda env list | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists. Activating..."
else
    echo "Environment '$ENV_NAME' does not exist. Creating and initializing..."
    # Create the environment from the YAML file
    conda env create -f slurm_environment.yml
fi
source activate base
conda init zsh
conda activate "$ENV_NAME"

configfile=$1

snakemake \
--profile simple/ \
--configfile $configfile 








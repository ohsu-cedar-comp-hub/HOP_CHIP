#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --partition batch
#SBATCH --time=8:00:00

ENV_NAME="HOP"

if conda env list | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists. Activating..."
else
    echo "Environment '$ENV_NAME' does not exist. Creating and initializing..."
    # Create the environment from the YAML file
    conda env create -f HOP_env.yaml
fi

eval "$(conda shell.bash hook)"
conda init
conda activate /home/exacloud/gscratch/CEDAR/chaoe/miniconda3/envs/HOP


# Default values
config=""
dir=""
profile=""
# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c) config="$2"; shift ;;
        -d) dir="$2"; shift ;;
        -p) profile="$2"; shift ;;

        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$config" ]; then
    echo "Error: -configfile is required"
    exit 1
fi

if [ -z "$dir" ] && [ -z "$profile" ]; then
    echo "No directory or profile provided to pull from. No S3 sync."
    
else 
    echo "Data directory and aws profile is provided. Syncing from S3..."
    aws --profile="$profile" --endpoint=https://rgw.ohsu.edu s3 sync s3://cedar-user-archive/chip_hop_test "$dir"
fi


snakemake --profile simple/ --configfile $config --cores 4  







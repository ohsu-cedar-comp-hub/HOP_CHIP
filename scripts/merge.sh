#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

in=""

cp $in $SCRATCH_PATH 

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

# 
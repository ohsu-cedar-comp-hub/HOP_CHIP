#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --gres=disk:1024 
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=2


# srun mkdir-scratch.sh
# SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --r1=*) r1="${1#*=}";;
        --r2=*) r2="${1#*=}";;
        --outdir=*) out="${1#*=}";;
        --threads=*) t="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$r1" || -z "$r2" || -z "$out" || -z "$t" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$out" ]; then
    mkdir -p "$out"
fi

#cp $r1 $SCRATCH_PATH
#cp $r2 $SCRATCH_PATH

echo "Running FASTQC:"
echo "R1:" $r1
echo "R2:" $r2
echo "Output Dir:" $out


fastqc $r1 \
    $r2 \
    -t $t \
    -o $out 

# srun rmdir-scratch.sh


#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

script=""
in=""
out=""


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --script=*) script="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$script" || -z "$out" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

lib=$(basename "$in" | cut -d. -f1)

cp $script $SCRATCH_PATH
cp $in $SCRATCH_PATH

annov_dir=$(dirname "$out") 

if [ ! -d "$annov_dir" ]; then
    mkdir -p "$annov_dir"
fi


echo "Formatting to fit Annovar:" 
echo "Using Script:" $SCRATCH_PATH/$(basename "$script")
echo "Input:" $SCRATCH_PATH/$(basename "$in")
echo "Output:" $out

python -u $SCRATCH_PATH/$(basename "$script")  --infile $SCRATCH_PATH/$(basename "$in") 

mv $SCRATCH_PATH/*.freq.annovar_file $out

srun rmdir-scratch.sh


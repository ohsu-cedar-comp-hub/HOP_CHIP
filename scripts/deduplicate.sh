#!/usr/bin/env bash

#SBATCH --time=20:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch
#SBATCH --mem=30G


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --outdir=*) out="${1#*=}";;
        --lib=*) lib="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$out" || -z "$lib"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$out" ]; then
    mkdir -p "$out"
fi


cp $in $SCRATCH_PATH

echo "Deduplicating.."
echo "BAM:" $SCRATCH_PATH/$(basename "$in")
echo "Writing output to:" $out

filename=$out/${lib}.sort.bam
metrics=$out/${lib}.dedup_metrics.txt

picard MarkDuplicates \
    I=$in \
    O=$filename \
    M=$metrics \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

mv $out/${lib}.sort.bai $out/${lib}.sort.bam.bai 

srun rmdir-scratch.sh

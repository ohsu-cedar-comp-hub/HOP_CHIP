#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --partition=batch
#SBATCH --gres=disk:1024 


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --outdir=*) out="${1#*=}";;
        --lib=*) lib="${1#*=}";;
        --rg=*) rg="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$out" || -z "$lib" || -z "$rg" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$out" ]; then
    mkdir -p "$out"
fi

checked_file=$out/${lib}.rg_exists.txt

if [ "$rg" = "true" ]; then
    cp $in $SCRATCH_PATH

    echo "Adding RG using Info From: " "$lib" 

    IFS="_" read -r exp sample libid flowcell lane <<< "$lib"
    lanenum=${lane#L}

    RGID="${sample}.${flowcell}.${lanenum}"
    RGLB="${libid}"
    RGPL="ILLUMINA"
    RGSM="${sample}"
    RGPU="${flowcell}.${lanenum}"

    samtools addreplacerg \
    -r "ID:$RGID" \
    -r "LB:$RGLB" \
    -r "PL:$RGPL" \
    -r "SM:$RGSM" \
    -r "PU:$RGPU" \
    -o $SCRATCH_PATH/${lib}.tmp.bam \
    $SCRATCH_PATH/$(basename "$in")

    mv "$SCRATCH_PATH/${lib}.tmp.bam" $in

    touch $checked_file


else
    echo "RG exists already! No changes made" 
    touch $checked_file

fi


srun rmdir-scratch.sh
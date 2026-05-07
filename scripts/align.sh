#!/usr/bin/env bash

#SBATCH --time=12:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=10


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --r1=*) r1="${1#*=}";;
        --r2=*) r2="${1#*=}";;
        --genome=*) genome="${1#*=}";;
        --outdir=*) out="${1#*=}";;
        --lib=*) lib="${1#*=}";;
        --threads=*) t="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$r1" || -z "$r2" || -z "$genome" || -z "$out" || -z "$lib" || -z "$t" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$out" ]; then
    mkdir -p "$out"
fi

cp $r1 $SCRATCH_PATH
cp $r2 $SCRATCH_PATH

prefix="${genome%.fa}"
cp ${prefix}.* "$SCRATCH_PATH"
echo "Copied genome index files:"
ls -lh "$SCRATCH_PATH"


echo "Aligning to genome:" $SCRATCH_PATH/$(basename "$genome")
echo "Read 1:" $SCRATCH_PATH/$(basename "$r1")
echo "Read 2:" $SCRATCH_PATH/$(basename "$r2")
echo "Writing output to:" $out

filename=$out/${lib}.aligned.bam

bwa-mem2 mem -t $t $SCRATCH_PATH/$(basename "$genome") $SCRATCH_PATH/$(basename "$r1") $SCRATCH_PATH/$(basename "$r2") \
  | samtools sort -@ $t -o $filename


srun rmdir-scratch.sh


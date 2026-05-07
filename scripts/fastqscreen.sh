#!/usr/bin/env bash
#SBATCH --time=5:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=8

srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

r1=""
r2=""
outdir=""
conf=""
threads=""
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --r1=*) r1="${1#*=}";;
        --r2=*) r2="${1#*=}";;
        --conf=*) conf="${1#*=}";;
        --outdir=*) outdir="${1#*=}";;
        --threads=*) t="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$r1" || -z "$r2" || -z "$conf" || -z "$outdir" || -z "$t" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

cp $r1 $SCRATCH_PATH 
cp $r2 $SCRATCH_PATH 
cp $conf $SCRATCH_PATH 


if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

cd $SCRATCH_PATH

fastq_screen \
  --conf $SCRATCH_PATH/$(basename "$conf") \
  --outdir $outdir \
  --threads $t  \
  $SCRATCH_PATH/$(basename "$r1") \
  $SCRATCH_PATH/$(basename "$r2")





srun rmdir-scratch.sh

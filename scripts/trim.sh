#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=2


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --r1=*) r1="${1#*=}";;
        --r2=*) r2="${1#*=}";;
        --r1_adaptor=*) r1_adaptor="${1#*=}";;
        --r2_adaptor=*) r2_adaptor="${1#*=}";;
        --outdir=*) out="${1#*=}";;
        --lib=*) lib="${1#*=}";;
        --threads=*) t="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$r1" || -z "$r2" || -z "$r1_adaptor" || -z "$r2_adaptor" ||  -z "$out" ||  -z "$lib" ||  -z "$t" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$out" ]; then
    mkdir -p "$out"
fi

cp $r1 $SCRATCH_PATH
cp $r2 $SCRATCH_PATH

echo "Trimming..:"
echo "R1:" $SCRATCH_PATH/$(basename "$r1")
echo "R2:" $SCRATCH_PATH/$(basename "$r2")
echo "Output Dir:" $out

cd $SCRATCH_PATH

cutadapt \
    -j $t \
    -a $r1_adaptor \
    -A $r2_adaptor \
    -o $SCRATCH_PATH/${lib}_trimmed_1.fq.gz \
    -p $SCRATCH_PATH/${lib}_trimmed_2.fq.gz \
    $SCRATCH_PATH/$(basename "$r1") $SCRATCH_PATH/$(basename "$r2")

mv $SCRATCH_PATH/${lib}_trimmed_1.fq.gz $out/${lib}_trimmed_1.fq.gz
mv $SCRATCH_PATH/${lib}_trimmed_2.fq.gz $out/${lib}_trimmed_2.fq.gz 

srun rmdir-scratch.sh
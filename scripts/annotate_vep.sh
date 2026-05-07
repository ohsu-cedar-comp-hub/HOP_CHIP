#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=6


# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --vep_dir=*) vep_dir="${1#*=}";;
        --t=*) t="${1#*=}";;
        --build=*) build="${1#*=}";;
        --outdir=*) outdir="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$vep_dir" || -z "$t" || -z "$build" || -z "$outdir" || -z "$out"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi



vep \
    -i "$in" \
    -o "$out" \
    --vcf \
    --cache \
    --dir_cache "$vep_dir" \
    --assembly "$build" \
    --everything \
    --check_existing \
    --compress_output  bgzip \
    --fork "$t" \
    --force_overwrite

echo "Finished annotation." 

bcftools index "$out"
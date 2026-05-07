#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1


# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --min_var_q=*) min_var_q="${1#*=}";;
        --min_depth=*) min_depth="${1#*=}";;
        --min_gq=*) min_gq="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$min_var_q" || -z "$in" || -z "$min_depth" || -z "$min_gq" || -z "$out"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi


lib=$(basename "$in" .raw.vcf.gz)

echo "Filtering via bcftools:"
echo "Input BAM:" "$in"
echo "Min Variant Quality Score:" "$min_var_q" 
echo "Min Depth:" "$min_depth" 
echo "Min Genotype Quality:" "$min_gq" 
echo "Output:" "$out"

# filter calls by depth and quality (bcftools)
bcftools filter \
    -e "QUAL<${min_var_q} || INFO/DP<${min_depth} || GQ<${min_gq}" \
    "$in" \
    -Oz -o "$out"


bcftools index "$out"

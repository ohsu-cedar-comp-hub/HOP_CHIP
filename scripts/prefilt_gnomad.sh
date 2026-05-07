#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1


# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --a=*) annotated_gnomad="${1#*=}";;
        --af=*) max_af="${1#*=}";;
        --in=*) in="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$annotated_gnomad" || -z "$in" || -z "$max_af" || -z "$out"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

lib=$(basename "$in")
lib=${lib%%_Var*.vcf.gz}

outdir=$(dirname "$out")

echo $lib

# prefilter germline variants: Annotate with gnomad and filter by population allelic fraction
# Note: VEP annotation below will also provide this annotation. Use this step only if you
# want to reduce the total number of calls to annotate with VEP by first filtering out
# variants that are likely germline
bcftools annotate \
    -a "$annotated_gnomad" \
    -c INFO/AF \
    "$in" \
    -Oz -o "$outdir/${lib}.AF_temp.vcf.gz"

echo "Finished annotating: $outdir/${lib}.AF_temp.vcf.gz"

bcftools index "$outdir/${lib}.AF_temp.vcf.gz"

echo "Prefiltering by Max AF: $max_af" 

bcftools filter \
    -e "INFO/AF>${max_af}" \
    "$outdir/${lib}.AF_temp.vcf.gz" \
    -Oz -o "$out"


bcftools index "$out"
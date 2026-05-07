#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1


# define args

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --regions=*) regions="${1#*=}";;
        --ref=*) ref="${1#*=}";;
        --in=*) in="${1#*=}";;
        --min_map_q=*) min_map_q="${1#*=}";;
        --min_base_q=*) min_base_q="${1#*=}";;
        --outdir=*) outdir="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$regions" || -z "$ref" || -z "$min_map_q" || -z "$in" || -z "$min_base_q" || -z "$outdir"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

lib=$(basename "$in" .sort.bam)

echo "Calling variants via bcftools:"
echo "Input BAM:" "$in"
echo "Reference:" "$ref"
echo "Regions:" "$regions"
echo "Min Map Quality:" "$min_map_q" 
echo "Min Base Quality:" "$min_base_q" 


# call SNV/indel variants (bcftools)
bcftools mpileup \
    -f "$ref" \
    -a AD,DP,SP \
    -q "$min_map_q" -Q "$min_base_q" \
    -T "$regions" \
    -Ou "$in" \
| bcftools call \
    -mv \
    -f GQ,GP \
    -Ou \
| bcftools norm \
    -f "$ref" \
    -m -any \
    -Oz -o "$outdir/${lib}.raw.vcf.gz"


bcftools index "$outdir/${lib}.raw.vcf.gz"

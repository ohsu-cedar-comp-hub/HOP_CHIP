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
        --a=*) annotated_gnomad="${1#*=}";;
        --num2chr=*) num2chr="${1#*=}";;
        --chr2num=*) chr2num="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$annotated_gnomad" || -z "$out"  || -z "$num2chr" || -z "$chr2num"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

# check if chrom names are aligned 
# if not change query vcf to match 

detect_style () {
    bcftools query -f '%CHROM\n' "$1" \
    | head -n 1000 \
    | grep -m 1 '^chr' >/dev/null && echo "chr" || echo "num"
}

style_in=$(detect_style "$in")
style_a=$(detect_style "$annotated_gnomad")

echo "CHROM in $in: $style_in"
echo "CHROM in $annotated_gnomad:  $style_a"


if [ "$style_in" != "$style_a" ]; then
    echo "CHROM mismatch detected — harmonizing $in to match $a"

    if [ "$style_a" = "chr" ]; then
            bcftools annotate \
            --rename-chrs "$num2chr" \
            -Oz -o "$out" \
            "$in"
    else
            bcftools annotate \
            --rename-chrs "$chr2num" \
            -Oz -o "$out" \
            "$in"
    fi 


else
    echo "CHROM formats already match, no conversion needed"
    cp $in $out 
fi









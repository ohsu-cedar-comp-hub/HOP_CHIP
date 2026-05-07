#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

min_depth=0
min_vaf=0
in=""
out=""
metrics=""
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --min_depth=*) min_depth="${1#*=}";;
        --min_vaf=*) min_vaf="${1#*=}";;
        --in=*) in="${1#*=}";;
        --metrics=*) metrics="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$metrics" || -z "$out"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

cp $in $SCRATCH_PATH 
cp $metrics $SCRATCH_PATH

readcounts_dir=$(dirname "$out") 
metrics_dir=$(dirname "$metrics") 

if [ ! -d "$readcounts_dir" ]; then
    mkdir -p "$readcounts_dir"
fi

if [ ! -d "$metrics_dir" ]; then
    mkdir -p "$metrics_dir"
fi


echo "Performing Filtering:"
echo "Input:" $SCRATCH_PATH/$(basename "$in")
echo "Min Depth:" $min_depth
echo "Min VAF:" $min_vaf
echo "Output:" $out 
echo "Metrics:" $metrics

# depth = 7th column 
# vaf = 6th column 
awk -v min_depth="$min_depth" -v min_vaf="$min_vaf" 'BEGIN { OFS = "\t" } $7 >= min_depth && $6 >= min_vaf { print $0 }' $SCRATCH_PATH/$(basename "$in") > $SCRATCH_PATH/$(basename "$out")

num=$(wc -l < "$SCRATCH_PATH/$(basename "$out")")

tmp=$(mktemp)

if [ "$(wc -l < $SCRATCH_PATH/$(basename "$metrics"))" -ge 3 ]; then
  # replace line 3 if exists 
  awk -v n="$num" -v d="$(date +"%Y-%m-%d %H:%M:%S")" -v r="$RULE_NAME" 'NR==3{$0=n "\t" d "\t" r} 1' "$SCRATCH_PATH/$(basename "$metrics")" > "$tmp" \
        && mv "$tmp" "$SCRATCH_PATH/$(basename "$metrics")"
else
  echo -e "$num\t$(date +"%Y-%m-%d %H:%M:%S")\t$RULE_NAME" >>  $SCRATCH_PATH/$(basename "$metrics")
fi


mv $SCRATCH_PATH/$(basename "$out") "$readcounts_dir"
mv $SCRATCH_PATH/$(basename "$metrics") "$metrics_dir"

srun rmdir-scratch.sh
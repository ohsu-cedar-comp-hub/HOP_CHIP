#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

script=""
in=""
exonic_func=""
out=""
metrics=""


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --script=*) script="${1#*=}";;
        --annot=*) annot="${1#*=}";;
        --metrics=*) metrics="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$script" || -z "$annot" || -z "$metrics" || -z "$out" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

cp $script $SCRATCH_PATH
cp $in $SCRATCH_PATH
cp $annot $SCRATCH_PATH
cp $metrics $SCRATCH_PATH

readcounts_dir=$(dirname "$in") 
annov_dir=$(dirname "$annot") 
metrics_dir=$(dirname "$metrics")
annotate_dir=$(dirname "$out")


if [ ! -d "$readcounts_dir" ]; then
    mkdir -p "$readcounts_dir"
fi

if [ ! -d "$annov_dir" ]; then
    mkdir -p "$annov_dir"
fi

if [ ! -d "$metrics_dir" ]; then
    mkdir -p "$metrics_dir"
fi

if [ ! -d "$annotate_dir" ]; then
    mkdir -p "$annotate_dir"
fi

echo "Adding Annotations:" 
echo "Using Script:" $SCRATCH_PATH/$(basename "$script")
echo "Input:" $SCRATCH_PATH/$(basename "$in")
echo "Annotation Info:" $SCRATCH_PATH/$(basename "$annot")
echo "Output:" $out
echo "Metrics:" $metrics


python -u $SCRATCH_PATH/$(basename "$script") --infile $SCRATCH_PATH/$(basename "$in") --annov $SCRATCH_PATH/$(basename "$annot") --outdir $SCRATCH_PATH

num=$(wc -l < $SCRATCH_PATH/*.annotated_readcounts)

tet2=$(awk -F'\t' '$13 == "TET2"' $SCRATCH_PATH/*.annotated_readcounts | wc -l)

dnmt3a=$(awk -F'\t' '$13 == "DNMT3A"' $SCRATCH_PATH/*.annotated_readcounts | wc -l)

tmp=$(mktemp)

timestamp=$(date +"%Y-%m-%d %H:%M:%S")

awk -v n="$num" -v t="$tet2" -v d3a="$dnmt3a" -v d="$timestamp" -v r="$RULE_NAME" '
{
  if (NR == 4) print n "\t" d "\t" r
  else if (NR == 5) print t "\t" d "\t" r "\t" "tet2"
  else if (NR == 6) print d3a "\t" d "\t" r "\t" "dnmt3a"
  else print
}
END {
  # append missing lines if total lines < 6
  for (i=NR+1; i<=6; i++) {
    if (i == 4) print n "\t" d "\t" r
    else if (i == 5) print t "\t" d "\t" r "\t" "tet2"
    else if (i == 6) print d3a "\t" d "\t" r "\t" "dnmt3a"
  }
}
' $SCRATCH_PATH/$(basename "$metrics") > $tmp && mv $tmp $SCRATCH_PATH/$(basename "$metrics")


mv $SCRATCH_PATH/*.annotated_readcounts "$annotate_dir"
mv $SCRATCH_PATH/$(basename "$metrics") "$metrics_dir"


srun rmdir-scratch.sh


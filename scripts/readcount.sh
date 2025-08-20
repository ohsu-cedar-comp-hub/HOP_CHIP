#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

max_depth=""
regions=""
ref=""
in=""
out=""
metrics=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --max_depth=*) max_depth="${1#*=}";;
        --regions=*) regions="${1#*=}";;
        --ref=*) ref="${1#*=}";;
        --in=*) in="${1#*=}";;
        --metrics=*) metrics="${1#*=}";;
        --out=*) out="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$max_depth" || -z "$regions" || -z "$ref" || -z "$in" || -z "$metrics" || -z "$out" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi




cp $regions $SCRATCH_PATH
cp $ref $SCRATCH_PATH
cp $in $SCRATCH_PATH
cp $in".bai" $SCRATCH_PATH

echo "Running bam-readcount using:"
echo "Max Depth:" $max_depth
echo "List of Regions:" $SCRATCH_PATH/$(basename "$regions")
echo "Reference:" $SCRATCH_PATH/$(basename "$ref")
echo "Input:" $SCRATCH_PATH/$(basename "$in")
echo "Writing output to:" $out
echo "Creating and adding to Metrics File:" $metrics


bam-readcount -d $max_depth \
        -w 1 \
        -l $SCRATCH_PATH/$(basename "$regions") \
        -f $SCRATCH_PATH/$(basename "$ref") \
        $SCRATCH_PATH/$(basename "$in")  > $SCRATCH_PATH/output


lib=$(basename "$in" | cut -d. -f1)

echo "Adding Library name as column:" $lib

awk -v lib="$lib" 'BEGIN{ FS = OFS = "\t" } { print lib, $0 }' \
    $SCRATCH_PATH/output > $SCRATCH_PATH/$(basename "$out")

num=$(wc -l < "$SCRATCH_PATH/$(basename "$out")")



tmp=$(mktemp)

if [ -f "$metrics" ]; then
    cp $metrics $SCRATCH_PATH
    awk -v n="$num" -v d="$(date +"%Y-%m-%d %H:%M:%S")" -v r="$RULE_NAME" 'NR==1{$0=n "\t" d "\t" r } 1' "$SCRATCH_PATH/$(basename "$metrics")" > "$tmp" \
        && mv "$tmp" "$SCRATCH_PATH/$(basename "$metrics")"
else
    echo -e "$num\t$(date +"%Y-%m-%d %H:%M:%S")\t$RULE_NAME" > "$SCRATCH_PATH/$(basename "$metrics")"
fi


mv $SCRATCH_PATH/$(basename "$out") $(dirname "$out") 
mv $SCRATCH_PATH/$(basename "$metrics") $(dirname "$out") 

srun rmdir-scratch.sh
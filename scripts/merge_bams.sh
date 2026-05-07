#!/usr/bin/env bash
#SBATCH --time=36:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1       
#SBATCH --gres=disk:1024          
#SBATCH --cpus-per-task=8


SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_PATH"

# -------------------------
# Parse arguments
# -------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dir1=*) dir1="${1#*=}";;
        --dir2=*) dir2="${1#*=}";;
        --outdir=*) outdir="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

[[ -z "${dir1:-}" || -z "${dir2:-}" || -z "${outdir:-}" ]] && {
    echo "Missing required arguments"
    exit 1
}

mkdir -p "$outdir"

THREADS=${SLURM_CPUS_PER_TASK:-1}

echo $THREADS

# -------------------------
# Copy BAMs safely
# -------------------------
echo "Copying BAMs to scratch..."

bams=( "$dir1"/*.sort.bam "$dir2"/*.sort.bam )

if (( ${#bams[@]} == 0 )); then
    echo "No BAM files found"
    exit 1
fi

for f in "${bams[@]}"; do
    cp -- "$f" "$SCRATCH_PATH/"
done

scratch_bams=( "$SCRATCH_PATH"/*.sort.bam )

# -------------------------
# Group BAMs by prefix (FAST + SAFE)
# -------------------------
declare -A groups

for f in "${scratch_bams[@]}"; do
    base=$(basename "$f")
    prefix="${base%%_*}_${base#*_}"
    prefix="${prefix%%_*}_${prefix#*_}"   # ensures first two fields only

    # safer extraction (first two underscore fields)
    prefix=$(awk -F'_' '{print $1"_"$2}' <<< "$base")

    groups["$prefix"]+="$f"$'\n'
done

# -------------------------
# Process each group
# -------------------------
for prefix in "${!groups[@]}"; do
    echo "Processing: $prefix"

    IFS=$'\n' read -r -d '' -a files \
        <<< "${groups[$prefix]}"$'\0'

    outbam="${SCRATCH_PATH}/${prefix}_merged.bam"
    final="${outdir}/${prefix}_merged.bam"

    if (( ${#files[@]} > 1 )); then
        echo "Merging ${#files[@]} files"
        samtools merge -@ "$THREADS" -o "$outbam" "${files[@]}"
    else
        echo "Single file → copying"
        cp -- "${files[0]}" "$outbam"
    fi

    echo "Indexing"
    samtools index -@ "$THREADS" "$outbam"

    mv -- "$outbam" "$final"
    mv -- "${outbam}.bai" "${final}.bai"
done

echo "Done."

srun rmdir-scratch.sh
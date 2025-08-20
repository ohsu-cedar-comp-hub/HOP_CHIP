#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch

srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"

dir=""
metadata=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dir=*) dir="${1#*=}";;
        --metadata=*) metadata="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$dir" || -z "$metadata"  ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

cp $dir/*.metrics $SCRATCH_PATH
cp $metadata $SCRATCH_PATH


awk -v dir="$SCRATCH_PATH" '
BEGIN {
  FS = OFS = "\t"
}

NR == 1 {
  # Print header with new columns
  print $0, "NucsPos", "VariantsPostQC", "VariantsFiltered", "AnnotatedVariants", "TET2", "DNM3TA"
  next
}

{
  lib = $10
  sub(".*/", "", lib)
  sub(/\.bam$/, "", lib)

  file = dir "/" lib ".metrics"

  # by default all metrics NA
  for (i = 1; i <= 6; i++) m[i] = "NA"

  # look for corresponding metrics file and fill in metrics for that lib row
  if ((getline line < file) > 0) {
    split(line, fields, "\t")
    m[1] = fields[1]
    n = 1

    while (n < 6 && (getline line < file) > 0) {
      n++
      split(line, fields, "\t")
      m[n] = fields[1]
    }
    close(file)
  } else {
    # if file missing or empty, just keep m[] as "NA"
    close(file)
  }

  print $0, m[1], m[2], m[3], m[4], m[5], m[6]
}
' "$SCRATCH_PATH/$(basename "$metadata")" > "$SCRATCH_PATH/enhanced_metadata"

mv "$SCRATCH_PATH/enhanced_metadata" "$dir/"
mv $SCRATCH_PATH/enhanced_metadata $dir 

srun rmdir-scratch.sh











#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --gres=disk:1024 
#SBATCH --partition=batch


srun mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
# define args

in=""
annovar=""
db=""
build=""
outdir=""
addl=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --in=*) in="${1#*=}";;
        --annovar=*) annovar="${1#*=}";;
        --db=*) db="${1#*=}";;
        --addl=*) addl="${1#*=}";;
        --build=*) build="${1#*=}";;
        --outdir=*) outdir="${1#*=}";;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

if [[ -z "$in" || -z "$annovar" || -z "$db" || -z "$build" || -z "$outdir" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi

annovarinfo_dir=$(dirname "$annovar")

cp $in $SCRATCH_PATH
cp $annovarinfo_dir/*.pl $SCRATCH_PATH
cp -r $db $SCRATCH_PATH



if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi


default_protocol="refGene"
operations="g"


if [[ -n "$addl" ]]; then

    IFS=',' read -ra addl_array <<< "$addl"

    protocol="$default_protocol,$addl"
    
    # generate additional f's 
    for _ in "${addl_array[@]}"; do
        operations+=",f"

    done

else
    protocol="$default_protocol"


fi

echo "Running annovar using:"
echo "Build:" $build
echo "Input:" $SCRATCH_PATH/$(basename "$in")
echo "DB:" $SCRATCH_PATH/$(basename "$db")
echo "DB Specifics:" $protocol 
echo "With These Operations:" $operations
echo "OutputDir:" $outdir

lib=$(basename $in | cut -d'.' -f1)

$SCRATCH_PATH/$(basename "$annovar") \
  $SCRATCH_PATH/$(basename "$in") \
  $SCRATCH_PATH/$(basename "$db") \
  -buildver $build \
  -out $SCRATCH_PATH/$lib \
  -protocol $protocol \
  -operation $operations \
  -remove \
  -nastring . \
  -polish

mv $SCRATCH_PATH/*.txt $outdir

srun rmdir-scratch.sh



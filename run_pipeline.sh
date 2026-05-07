#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --ntasks=1 
#SBATCH --mem 4G
#SBATCH --partition batch
#SBATCH --time=36:00:00


# Default values
config=""
dir=""
profile=""
until=""
cores=6
env=""

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c) config="$2"; shift ;;
        -env) env="$2"; shift ;;
        -d) dir="$2"; shift ;;
        -s3_p) s3profile="$2"; shift ;;
        -clust_p) cluster_p="$2"; shift ;;
        -u) until="$2"; shift ;;
        -cores) cores="$2"; shift ;;
        -s) sfile="$2"; shift ;;
        -f) forcerun="$2"; shift ;;

        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$config" ]; then
    echo "Error: -configfile is required"
    exit 1
fi

if [ -z "$env" ]; then
    echo "Error: -env is required"
    exit 1
fi

if conda env list | awk '{print $1}' | grep -qx "$env"; then
    echo "Environment '$env' already exists. Activating..."

    eval "$(conda shell.bash hook)"
    conda activate "$env"

else
    echo "Environment '$env' doesn't exist. Please install first."
    exit 1
fi

if [ -z "$dir" ] && [ -z "$s3profile" ]; then
    echo "No directory or profile provided to pull from. No S3 sync."
    
else 
    echo "Data directory and aws profile is provided. Syncing from S3..."
    aws --profile="$profile" --endpoint=https://rgw.ohsu.edu s3 sync s3://cedar-user-archive/chip_hop_test "$dir"
fi

cmd="snakemake --profile $cluster_p --configfile $config --cores $cores "

if [ -n "$until" ]; then
    cmd="$cmd --until $until"
fi

if [ -n "$forcerun" ]; then
    cmd="$cmd --forcerun $forcerun"
fi


if [ -n "$sfile" ]; then
    cmd="$cmd -s $sfile"
fi


echo "$cmd" 

eval "$cmd"









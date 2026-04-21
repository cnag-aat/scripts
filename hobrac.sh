#!/bin/bash

# Recommended sbatch settings:
# --time=06:00:00
# --cpus-per-task=8
# --mem=128G
ASSEMBLY=$1
SPECIES=$2
TAXID=$3

# Set environment variable
export TAXONKIT_DB=/software/assembly/containers/hobrac/taxonkit_db

# Optional: avoid cache issues on shared systems
export SINGULARITY_CACHEDIR=/scratch_isilon/groups/assembly/$USER/singularity_cache

# Input/output paths (edit if needed)
SIF=/software/assembly/containers/hobrac/hobrac-tools_0.1.5.sif
WORKDIR=$PWD

echo "Starting HOBRAC job on $(date)"
echo "Running in $WORKDIR"

# Run command inside container
/usr/bin/singularity exec \
    --bind $WORKDIR:$WORKDIR \
    --pwd $WORKDIR \
    --bind /software/assembly/containers/hobrac/taxonkit_db:/taxonkit_db \
    --env TAXONKIT_DB=/taxonkit_db \
    $SIF \
    hobrac \
        -a $ASSEMBLY \
        -n "$SPECIES" \
        -t $TAXID \
        -o hobrac_out \
        -e local

echo "Finished on $(date)"

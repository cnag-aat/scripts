#!/usr/bin/env bash

# |-- ENA
# |-- analysis
# |-- annotation
# |-- assembly
# |   |-- curation
# |   |-- hic_qc
# |   `-- mitogenome
# `-- data
#     |-- cnag
#     |-- external
#     `-- reference

SPECIES=$1
CNAGDATA=$2
BASE=$PWD;
mkdir -p $BASE/$SPECIES
chmod 775 $BASE/$SPECIES
mkdir -p "$BASE/$SPECIES/assembly/asm1"
chmod 775 "$BASE/$SPECIES/assembly/asm1"
mkdir -p "$BASE/$SPECIES/analysis"
chmod 775 "$BASE/$SPECIES/analysis"
mkdir -p "$BASE/$SPECIES/annotation"
chmod 775 "$BASE/$SPECIES/annotation"
mkdir -p "$BASE/$SPECIES/assembly/curation"
chmod 775 "$BASE/$SPECIES/assembly/curation"
mkdir -p "$BASE/$SPECIES/assembly/mitogenome"
chmod 775 "$BASE/$SPECIES/assembly/mitogenome"
mkdir -p "$BASE/$SPECIES/assembly/cobiont"
chmod 775 "$BASE/$SPECIES/assembly/cobiont"
mkdir -p "$BASE/$SPECIES/assembly/hic_qc"
chmod 775 "$BASE/$SPECIES/assembly/hic_qc"
mkdir -p "$BASE/$SPECIES/data"
chmod 775 "$BASE/$SPECIES/data"
if [ $2 ]; then
    ln -s $CNAGDATA "$BASE/$SPECIES/data/cnag"
fi
mkdir -p "$BASE/$SPECIES/data/external"
chmod 775 "$BASE/$SPECIES/data/external"
mkdir -p "$BASE/$SPECIES/data/reference"
chmod 775 "$BASE/$SPECIES/data/reference"
chmod 775 "$BASE/$SPECIES/data"
mkdir -p "$BASE/$SPECIES/ENA/"
chmod 775 "$BASE/$SPECIES/ENA"
chmod -R +t "$BASE/$SPECIES"


#!/usr/bin/env bash


# |-- analysis
# |-- annotation
# |-- assembly
# |-- curation
# |-- data
# |   |-- dna
# |   |   |-- clr
# |   |   |-- hic
# |   |   |-- hifi
# |   |   |-- illumina
# |   |   `-- ont
# |   |-- external
# |   |-- reference
# |   `-- rna
# |       |-- illumina
# |       |-- isoseq
# |       `-- ont
# `-- preqc

SPECIES=$1
CNAGDATA=$2
BASE=$PWD;
mkdir -p $BASE/$SPECIES
chmod 775 $BASE/$SPECIES
mkdir -p "$BASE/$SPECIES/assembly"
chmod 775 "$BASE/$SPECIES/assembly"
mkdir -p "$BASE/$SPECIES/analysis"
chmod 775 "$BASE/$SPECIES/analysis"
mkdir -p "$BASE/$SPECIES/annotation"
chmod 775 "$BASE/$SPECIES/annotation"
mkdir -p "$BASE/$SPECIES/assembly/curation"
chmod 775 "$BASE/$SPECIES/assembly/curation"
mkdir -p "$BASE/$SPECIES/assembly/mitogenome"
chmod 775 "$BASE/$SPECIES/assembly/mitogenome"
mkdir -p "$BASE/$SPECIES/assembly/hic_qc"
chmod 775 "$BASE/$SPECIES/assembly/hic_qc"
mkdir -p "$BASE/$SPECIES/data/cnag"
chmod 775 "$BASE/$SPECIES/data/cnag"
for i in $CNAGDATA/*; do ln -s $i $BASE/$SPECIES/data/cnag/; done
mkdir -p "$BASE/$SPECIES/data/external"
chmod 775 "$BASE/$SPECIES/data/external"
mkdir -p "$BASE/$SPECIES/data/reference"
chmod 775 "$BASE/$SPECIES/data/reference"
chmod 775 "$BASE/$SPECIES/data"
mkdir -p "$BASE/$SPECIES/ENA/"
chmod 775 "$BASE/$SPECIES/ENA"
chmod -R +t "$BASE/$SPECIES"


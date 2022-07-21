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
BASE='/scratch/devel/talioto/denovo_assemblies'
mkdir -p $BASE/$SPECIES
chmod 775 $BASE/$SPECIES
mkdir -p "$BASE/$SPECIES/assembly"
chmod 775 "$BASE/$SPECIES/assembly"
mkdir -p "$BASE/$SPECIES/analysis"
chmod 775 "$BASE/$SPECIES/analysis"
mkdir -p "$BASE/$SPECIES/annotation"
chmod 775 "$BASE/$SPECIES/annotation"
mkdir -p "$BASE/$SPECIES/curation"
chmod 775 "$BASE/$SPECIES/curation"
mkdir -p "$BASE/$SPECIES/data/dna/ont"
chmod 775 "$BASE/$SPECIES/data/dna/ont"
mkdir -p "$BASE/$SPECIES/data/dna/illumina"
chmod 775 "$BASE/$SPECIES/data/dna/illumina"
mkdir -p "$BASE/$SPECIES/data/dna/hic"
chmod 775 "$BASE/$SPECIES/data/dna/hic"
mkdir -p "$BASE/$SPECIES/data/rna/ont"
chmod 775 "$BASE/$SPECIES/data/rna/ont"
mkdir -p "$BASE/$SPECIES/data/rna/illumina"
chmod 775 "$BASE/$SPECIES/data/rna/illumina"
mkdir -p "$BASE/$SPECIES/data/external"
chmod 775 "$BASE/$SPECIES/data/external"
mkdir -p "$BASE/$SPECIES/data/reference"
chmod 775 "$BASE/$SPECIES/data/reference"
chmod 775 "$BASE/$SPECIES/data"
chmod 775 "$BASE/$SPECIES/data/dna"
chmod 775 "$BASE/$SPECIES/data/rna"

mkdir -p "$BASE/$SPECIES/preqc"
chmod 775 "$BASE/$SPECIES/preqc"

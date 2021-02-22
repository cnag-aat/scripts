#!/bin/bash

READS=$1  #reads.fastq
BASE=$2   #SAMPLEBC_v5
THREADS=$3  #16
GSIZE=$4  #5.3m
HYBRID=$5 #hybrid assembly
##################
###### FLYE ######
##################
source /home/devel/talioto/miniconda3/etc/profile.d/conda.sh
conda deactivate
if [ ! -f $BASE.flye.fasta ]; then
    module purge && module load gcc/6.3.0 openssl/1.0.2q  python/3.7.5 flye
    export PATH=/apps/FLYE/2.6/bin/:$PATH
    flye -o $BASE.FLYE -t $THREADS -i 2 -g $GSIZE --plasmids --meta --nano-raw  $READS
fi
##################
##### RACON ######
##################
if [ ! -f $BASE.flye.racon1.racon2.fasta ]; then
    ln -s $BASE.FLYE/assembly.fasta $BASE.flye.fasta
    module purge && module load minimap2 gcc/4.9.3-gold racon
    /scratch/project/devel/aateam/bin/polish.pl -o $READS -r 2 -p 0 -t $THREADS -i $BASE.flye.fasta
fi
##################
##### MEDAKA #####
##################
if [ ! -f $BASE.flye.racon.medaka.fasta ]; then
    module purge && conda activate medaka
    medaka_consensus -t $THREADS -m r941_min_high_g344 -i $READS -d $BASE.flye.racon1.racon2.fasta -o $TMPDIR/$BASE.flye.racon1.racon2.medaka;  cp -r $TMPDIR/$BASE.flye.racon1.racon2.medaka .; rm -r $TMPDIR/$BASE.flye.racon1.racon2.medaka

    conda deactivate
    ln -s $BASE.flye.racon1.racon2.medaka/consensus.fasta $BASE.flye.racon.medaka.fasta
    /project/devel/aateam/bin/scaffolds2contigs.pl -i $BASE.flye.racon.medaka.fasta -name $BASE
fi

dnadiff -p $BASE.DNADIFF $HYBRID $BASE.scaffolds.fa

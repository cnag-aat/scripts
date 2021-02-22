#!/bin/bash

ASSEMBLY=$1
R9READS=$2  #reads.fastq
#R10READS=$2
PE1READS=$3
PE2READS=$4
BASE=$5   #SAMPLEBC_v5
THREADS=$6  #16
source /home/devel/talioto/miniconda3/etc/profile.d/conda.sh
conda deactivate
##################
##### RACON ######
##################
if [ ! -f $BASE.racon1.racon2.fasta ]; then
    ln -s $ASSEMBLY $BASE.fasta
    module purge && module load minimap2 gcc/4.9.3-gold racon
    /scratch/project/devel/aateam/bin/polish.pl -o $R9READS -rr 2 -pr 0 -t $THREADS -i $BASE.fasta
fi
##################
##### MEDAKA #####
##################
if [ ! -f $BASE.racon.medaka.fasta ]; then
    module purge && conda activate /home/devel/talioto/miniconda3/envs/medaka
    medaka_consensus -t $THREADS -m r941_min_high_g344 -i $R9READS -d $BASE.racon1.racon2.fasta -o $TMPDIR/$BASE.racon1.racon2.medaka;  cp -r $TMPDIR/$BASE.racon1.racon2.medaka .; rm -r $TMPDIR/$BASE.racon1.racon2.medaka

    conda deactivate
    ln -s $BASE.racon1.racon2.medaka/consensus.fasta $BASE.racon.medaka.fasta
    /project/devel/aateam/bin/scaffolds2contigs.pl -i $BASE.racon.medaka.fasta -name $BASE
fi

##################
##### PILON ######
##################
if [ ! -f $BASE.racon.medaka.pilon.fasta ]; then
    module purge && module load gcc/4.9.3-gold racon bwa minimap2 samtools java/1.8.0u31 PILON
    /scratch/project/devel/aateam/bin/polish.pl --pe1 $PE1READS --pe2 $PE2READS -pr 2 -t $THREADS -i $BASE.racon.medaka.fasta
    ln -s $BASE.racon.medaka.pilon1.pilon2.fasta $BASE.racon.medaka.pilon.fasta
    /project/devel/aateam/bin/scaffolds2contigs.pl -i $BASE.racon.medaka.pilon.fasta -name $BASE
fi

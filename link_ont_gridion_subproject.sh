#!/usr/bin/env bash

export SP=$1
echo "copying $SP..."
for i in `find /scratch/production/shared/gridion/Runs/*$SP*/ -maxdepth 0 -mindepth 0`;
do export b=`basename $i`;
   mkdir $b; cd $b;
   rsync -ar --exclude=*fast5* --exclude=*fastq* --exclude=*fail* $i . >> ../rsync.$SP.out 2>>../rsync.$SP.error;
   ln -s $i/*fastq.gz .;
   ln -s $i/*fast5.tar.gz .;
   cd ..;
done

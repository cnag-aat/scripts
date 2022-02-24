#!/usr/bin/env bash

export SP=$1
echo "copying $SP..."
for i in `find /scratch/project/production/promethion01/*$SP*/ -maxdepth 3 -mindepth 3`;
do rsync -ar --exclude=*fast5* --exclude=*fastq* --exclude=*fail* $i . > rsync.$SP.out 2>rsync.$SP.error;
   export b=`basename $i`;
   echo $b;
   ln -s $i/fast5_pass $b/;
   ln -s $i/fastq_pass $b/;
   zcat $i/fastq_pass/*fastq.gz | pigz > $b/$b.pass.fastq.gz;
done

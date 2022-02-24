#!/usr/bin/env bash

export SP=$1
echo "copying $SP..."
for i in `find /scratch/project/production/promethion01/Runs/*$SP*/ -maxdepth 0 -mindepth 0`;
do export b=`basename $i`;
    echo $b;
    mkdir -p $b;
    cd $b;
    rsync -ar --exclude=*fast5* --exclude=*err --exclude=*out --exclude=*backup* --exclude=*list --exclude=*sbatch* $i . > rsync.$SP.out 2>rsync.$SP.error;
    ln -s /scratch/project/production/promethion01/Runs/$b/*fast5.tar.gz;
    cd ..;
done

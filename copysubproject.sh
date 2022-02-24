#!/usr/bin/env bash

export SP=$1
echo "copying $SP..."
for i in `find /scratch/project/production/promethion01/*$SP*/ -maxdepth 3 -mindepth 3`; do rsync -ar --exclude=*fast5* --exclude=*fail* $i . > rsync.$SP.out 2>rsync.$SP.error; done

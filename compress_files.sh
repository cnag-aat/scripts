#!/bin/bash

target=$1
echo "compressing $target files"
find . -user talioto -type f -name  "*.$target" -exec echo pigz {} \; | head -n 10000 > $target.$$.args 
submit_arrayjob.pl -a $target.$$.args -name compress_$target.$$ -t 1
#rm $target.$$.args

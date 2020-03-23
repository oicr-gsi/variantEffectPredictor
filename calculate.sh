#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
find . -name "*.vcf.gz" -exec md5sum {} +
count=`ls -1 *.maf.gz | wc -l`
if [ $count !=0 ]; then
   zcat *.maf.gz | cut -f1 | sort | uniq | wc -l
   zcat *.maf.gz | cut -f16 | sort | uniq | wc -l
fi

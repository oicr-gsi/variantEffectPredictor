#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
ls | sort

find . -name '*.vep.vcf.gz' | xargs zcat | grep -v ^# |  cut -f8 --complement | md5sum | sort
find . -name '*.targeted.vcf.gz' | xargs -i -n 1 --no-run-if-empty md5sum {} | sort
for i in $(find . -name '*.maf.gz'); do if [ -f $i ]; then zcat $i | cut -f46,51 --complement | md5sum; fi; done
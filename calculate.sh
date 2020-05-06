#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
ls | sort

find -name *.vcf.gz | xargs zcat | grep -v ^# | md5sum | sort
find -name *.maf.gz | md5sum | sort

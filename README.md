# variantEffectPredictor

Variant Effect Predictor Workflow version 2.1

## Overview

## Dependencies

* [bedtools 2.27](https://github.com/arq5x/bedtools)
* [tabix 0.2.6](https://github.com/samtools/tabix)
* [vep 92.0](https://github.com/Ensembl/ensembl-vep)
* [vcftools 0.1.16](https://vcftools.github.io/index.html)


## Usage

### Cromwell
```
java -jar cromwell.jar run variantEffectPredictor.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`vcfFile`|File|Input VCF file
`vcfIndex`|File|Input VCF index file
`toMAF`|Boolean|If true, generate the MAF file
`onlyTumor`|Boolean|If true, run tumor only mode
`vep.ncbiBuild`|String|The assembly version
`vep.vepCacheDir`|String|Directory of cache files
`vep.referenceFasta`|String|Reference fasta file
`vep.modules`|String|Required environment modules
`vcf2maf.modules`|String|Required environment modules
`vcf2maf.referenceFasta`|String|Reference fasta file
`vcf2maf.ncbiBuild`|String|The assembly version
`vcf2maf.vepPath`|String|Path to vep script
`vcf2maf.vepCacheDir`|String|Directory of vep cache files
`vcf2maf.vcfFilter`|String|Filter for the vep module that is used in vcf2maf


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBed`|String?|None|Target bed file


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBedTask.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`targetBedTask.modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`targetBedTask.jobMemory`|Int|32|Memory allocated for this job (GB)
`targetBedTask.threads`|Int|4|Requested CPU threads
`targetBedTask.timeout`|Int|6|Hours before task timeout
`getSampleNames.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`getSampleNames.modules`|String|"vcftools/0.1.16"|Required environment modules
`getSampleNames.jobMemory`|Int|32|Memory allocated for this job (GB)
`getSampleNames.threads`|Int|4|Requested CPU threads
`getSampleNames.timeout`|Int|6|Hours before task timeout
`chromosomeArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`chromosomeArray.threads`|Int|4|Requested CPU threads.
`chromosomeArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`subsetVcf.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`subsetVcf.modules`|String|"bcftools/1.9"|Required environment modules
`subsetVcf.jobMemory`|Int|32|Memory allocated to job (in GB).
`subsetVcf.threads`|Int|4|Requested CPU threads.
`subsetVcf.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`vep.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`vep.addParam`|String?|None|Additional vep parameters
`vep.species`|String|"homo_sapiens"|Species name
`vep.jobMemory`|Int|32|Memory allocated for this job (GB)
`vep.threads`|Int|4|Requested CPU threads
`vep.timeout`|Int|16|Hours before task timeout
`tumorOnlyAlign.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorOnlyAlign.modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`tumorOnlyAlign.jobMemory`|Int|32|Memory allocated for this job (GB)
`tumorOnlyAlign.threads`|Int|4|Requested CPU threads
`tumorOnlyAlign.timeout`|Int|6|Hours before task timeout
`tumorOnlyAlign.updateTagValue`|Boolean|false|If true, update tag values in vcf header for CC workflow
`vcf2maf.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`vcf2maf.species`|String|"homo_sapiens"|Species name
`vcf2maf.retainInfoProvided`|Boolean|false|Comma-delimited names of INFO fields to retain as extra columns in MAF
`vcf2maf.maxfilterAC`|Int|10|The maximum AC filter
`vcf2maf.minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`vcf2maf.bufferSize`|Int|200|The buffer size
`vcf2maf.jobMemory`|Int|32|Memory allocated for this job (GB)
`vcf2maf.threads`|Int|4|Requested CPU threads
`vcf2maf.timeout`|Int|48|Hours before task timeout
`mergeMafs.modules`|String|"tabix/0.2.6"|Required environment modules
`mergeMafs.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeMafs.threads`|Int|4|Requested CPU threads.
`mergeMafs.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mergeVcfs.modules`|String|"gatk/4.1.7.0"|Required environment modules.
`mergeVcfs.extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`mergeVcfs.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeVcfs.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeVcfs.threads`|Int|4|Requested CPU threads.
`mergeVcfs.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description
---|---|---
`outputVcf`|File|Annotated vcf output file from vep
`outputTbi`|File|Index of the annotated vcf output file from vep
`outputMaf`|File?|Maf output file from vcf2maf(if toMAF is true)
`outputTargetVcf`|File?|Vcf on target for the input vcf (if targetBed is given), non annotated
`outputTargetTbi`|File?|Index of the vcf on target for the input vcf (if targetBed is given), non annotated


## Commands
 This section lists command(s) run by variantEffectPredictor workflow
 
 * Running variantEffectPredictor
 
 
 A workflow for annotating the SNV and INDEL mutation calls in VCF format, and generating a MAF file of annotated calls.
 ### (needs description)
 ```
     set -euo pipefail
 
     bedtools intersect -header -u \
                        -a ~{vcfFile} \
                        -b ~{targetBed} \
                        > ~{basename}.targeted.vcf
     
     bgzip -c ~{basename}.targeted.vcf > ~{basename}.targeted.vcf.gz
            
     tabix -p vcf ~{basename}.targeted.vcf.gz
 ```
 ### (needs description)
 ```
     zcat ~{vcfFile} | grep -v ^# | cut -f 1 | uniq
 ```
 ### (needs description)
 ```
     set -euo pipefail
 
     bcftools view -r ~{regions} ~{vcfFile} | bgzip -c > ~{basename}.vcf.gz 
 ```
 ### (needs description)
 
 ```
     set -euo pipefail
 
     if [ "~{species}" = "homo_sapiens" ]; then
       human_only_command_line="--polyphen b --af --af_1kg --af_esp --af_gnomad"
     else
       human_only_command_line=""
     fi
 
     vep --offline --dir ~{vepCacheDir} -i ~{vcfFile} --fasta ~{referenceFasta} --species ~{species} \
           --assembly ~{ncbiBuild} -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip ~{addParam} \
           --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype \
           --canonical --protein --biotype --uniprot --tsl --variant_class --check_existing --total_length \
           --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele \
           --pick_order canonical,tsl,biotype,rank,ccds,length  \
           $human_only_command_line \
           --pubmed --fork 4 --regulatory
 
 ```
 ### (needs description)
 ```
     set -euo pipefail
 
     vcf-query -l  "~{vcfFile}" > sample_headers_all
     cat sample_headers_all | grep -v "GATK" | tr "\n" "," > sample_names_all
     if [[ `cat sample_names_all | tr "," "\n" | wc -l` == 2 ]]; then
       for item in `cat sample_names_all | tr "," "\n"`; do if [[ $item == "NORMAL" || $item == *_R_* || $item == *_R || $item == *BC*  || $item == "unmatched" ]]; then NORM=$item; else TUMR=$item; fi; done
     else TUMR=`cat sample_names_all | tr -d ","`; NORM="unmatched"; fi
 
     echo $TUMR > names.txt
     echo $NORM >> names.txt
 
 ```
 ### (needs description)
 ```
     set -euo pipefail
 
     if ~{updateTagValue} ; then
         zcat ~{vcfFile} | sed s/Number\=A/Number\=./ | sed s/Number\=R/Number\=./ > "~{basename}_temporary.vcf"
         cat ~{basename}_temporary.vcf | sed 's/QSS\,Number\=A/QSS\,Number\=\./' | sed 's/AS_FilterStatus\,Number\=A/AS_FilterStatus\,Number\=\./' | bgzip -c > "~{basename}_input.vcf.gz"
     else
         zcat ~{vcfFile} | sed 's/QSS\,Number\=A/QSS\,Number\=\./' | sed 's/AS_FilterStatus\,Number\=A/AS_FilterStatus\,Number\=\./' | bgzip -c > "~{basename}_input.vcf.gz"
     fi
     
     tabix -p vcf "~{basename}_input.vcf.gz"
 
     cat ~{tumorNormalNames} > "~{basename}_header"
     bcftools merge "~{basename}_input.vcf.gz" "~{basename}_input.vcf.gz" --force-samples > "~{basename}.temp_tumor.vcf"
     bcftools reheader -s "~{basename}_header" "~{basename}.temp_tumor.vcf" > "~{basename}.unmatched.vcf"
     bgzip -c "~{basename}.unmatched.vcf" > "~{basename}.unmatched.vcf.gz"
     tabix -p vcf "~{basename}.unmatched.vcf.gz"
 ```
 ### (needs description)
 ```
     set -euo pipefail
 
     TUMR=$(sed -n 1p ~{tumorNormalNames} )
     NORM=$(sed -n 2p ~{tumorNormalNames} )
 
     bgzip -c -d ~{vcfFile} > ~{basename}
 
     if ~{retainInfoProvided} ; then
 
         vcf2maf --ref-fasta ~{referenceFasta} --species ~{species} --ncbi-build ~{ncbiBuild} \
                 --input-vcf ~{basename} --output-maf ~{basename}.maf \
                 --tumor-id $TUMR --normal-id $NORM --vcf-tumor-id $TUMR --vcf-normal-id $NORM \
                 --filter-vcf ~{vcfFilter} --vep-path ~{vepPath} --vep-data ~{vepCacheDir} \
                 --max-filter-ac ~{maxfilterAC} --min-hom-vaf ~{minHomVaf} --buffer-size ~{bufferSize} --retain-info MBQ,MMQ,TLOD,set
     else
         vcf2maf --ref-fasta ~{referenceFasta} --species ~{species} --ncbi-build ~{ncbiBuild} \
                 --input-vcf ~{basename} --output-maf ~{basename}.maf \
                 --tumor-id $TUMR --normal-id $NORM --vcf-tumor-id $TUMR --vcf-normal-id $NORM \
                 --filter-vcf ~{vcfFilter} --vep-path ~{vepPath} --vep-data ~{vepCacheDir} \
                 --max-filter-ac ~{maxfilterAC} --min-hom-vaf ~{minHomVaf} --buffer-size ~{bufferSize}
     fi
 ```
 ### (needs description)
 ```
     set -euo pipefail
 
     head -n 2 ~{mafs[0]} > ~{basename}
     cat ~{sep=" " mafs} | grep -v ^# | grep -v "Hugo_Symbol" >> ~{basename}
     bgzip -c ~{basename} > ~{basename}.maf.gz
 
 ```
 ### (needs description)
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeVcfs \
     -I ~{sep=" -I " vcfs} ~{extraArgs} \
     -O ~{basename}.vcf.gz
 ```
  
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

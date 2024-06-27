# variantEffectPredictor

Variant Effect Predictor Workflow version 2.2

## Overview

## Dependencies

* [bedtools 2.27](https://github.com/arq5x/bedtools)
* [tabix 0.2.6](https://github.com/samtools/tabix)
* [vep 105.0](https://github.com/Ensembl/ensembl-vep)
* [vcf2maf 1.6.21b](https://github.com/mskcc/vcf2maf/commit/5ed414428046e71833f454d4b64da6c30362a89b)
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
`tumorName`|String|Name of the tumor sample
`toMAF`|Boolean|If true, generate the MAF file
`onlyTumor`|Boolean|If true, run tumor only mode
`reference`|String|reference genome for input sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBed`|String?|None|Target bed file
`normalName`|String?|None|Name of the normal sample


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBedTask.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`targetBedTask.modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`targetBedTask.jobMemory`|Int|12|Memory allocated for this job (GB)
`targetBedTask.timeout`|Int|6|Hours before task timeout
`getSampleNames.jobMemory`|Int|1|Memory allocated for this job (GB)
`getSampleNames.timeout`|Int|1|Hours before task timeout
`chromosomeArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`chromosomeArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`getChrCoefficient.memory`|Int|1|Memory allocated for this job
`getChrCoefficient.timeout`|Int|1|Hours before task timeout
`subsetVcf.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`subsetVcf.modules`|String|"bcftools/1.9"|Required environment modules
`subsetVcf.jobMemory`|Int|12|Memory allocated to job (in GB).
`subsetVcf.minMemory`|Int|4|Minimum RAM (in GB) allocated to the task
`subsetVcf.timeout`|Int|2|Maximum amount of time (in hours) the task can run for.
`vep.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`vep.addParam`|String?|None|Additional vep parameters
`vep.vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`vep.jobMemory`|Int|12|Memory allocated for this job (GB)
`vep.minMemory`|Int|4|Minimum RAM (in GB) allocated to the task
`vep.threads`|Int|4|Requested CPU threads
`vep.timeout`|Int|16|Hours before task timeout
`tumorOnlyAlign.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorOnlyAlign.modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`tumorOnlyAlign.jobMemory`|Int|12|Memory allocated for this job (GB)
`tumorOnlyAlign.minMemory`|Int|4|Minimum RAM (in GB) allocated to the task
`tumorOnlyAlign.timeout`|Int|6|Hours before task timeout
`tumorOnlyAlign.updateTagValue`|Boolean|false|If true, update tag values in vcf header for CC workflow
`vcf2maf.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`vcf2maf.retainInfoProvided`|Boolean|false|Comma-delimited names of INFO fields to retain as extra columns in MAF
`vcf2maf.vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`vcf2maf.minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`vcf2maf.bufferSize`|Int|200|The buffer size
`vcf2maf.jobMemory`|Int|12|Memory allocated for this job (GB)
`vcf2maf.minMemory`|Int|4|A minimum amount of memory allocated to the task, overrides the scaled RAM setting
`vcf2maf.threads`|Int|4|Requested CPU threads
`vcf2maf.timeout`|Int|18|Hours before task timeout
`mergeMafs.modules`|String|"tabix/0.2.6"|Required environment modules
`mergeMafs.jobMemory`|Int|8|Memory allocated to job (in GB).
`mergeMafs.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mergeVcfs.modules`|String|"gatk/4.1.7.0"|Required environment modules.
`mergeVcfs.extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`mergeVcfs.jobMemory`|Int|8|Memory allocated to job (in GB).
`mergeVcfs.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeVcfs.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description | Labels
---|---|---|---
`outputVcf`|File|Annotated vcf output file from vep|vidarr_label: outputVcf
`outputTbi`|File|Index of the annotated vcf output file from vep|vidarr_label: outputTbi
`outputMaf`|File?|Maf output file from vcf2maf(if toMAF is true)|vidarr_label: outputMaf
`outputTargetVcf`|File?|Vcf on target for the input vcf (if targetBed is given), non annotated|vidarr_label: outputTargetVcf
`outputTargetTbi`|File?|Index of the vcf on target for the input vcf (if targetBed is given), non annotated|vidarr_label: outputTargetTbi


## Commands
This section lists command(s) run by variantEffectPredictor workflow
 
* Running variantEffectPredictor
 
 
### Derive a chromosome-specific scaling coefficient
 
```
     CHR_LEN=$(zcat ~{inputVcf} | head -n 800 | grep contig | grep -w ~{chromosome} | grep -v _ | sed -r 's/.*length=([[:digit:]]+)./\1/')
     LARGEST=$(zcat ~{inputVcf} | head -n 800 | grep contig | grep -v _ | sed -r 's/.*length=([[:digit:]]+)./\1/' | sort -n | tail -n 1)
     echo | awk -v chr_len=$CHR_LEN -v largest_chr=$LARGEST '{print int((chr_len/largest_chr + 0.1) * 10)/10}'
```
 
### targetBedTask - extract variants overlapping intervals from a .bed file
 
```
     set -euo pipefail
 
     bedtools intersect -header -u \
                        -a ~{vcfFile} \
                        -b ~{targetBed} \
                        > ~{basename}.targeted.vcf
 
     bgzip -c ~{basename}.targeted.vcf > ~{basename}.targeted.vcf.gz
 
     tabix -p vcf ~{basename}.targeted.vcf.gz
```
 
### Retrieve all chromosomes from vcf file
 
This ensures that we split only by those chromosomes for which we have records
 
```
     zgrep -v ^# ~{vcfFile} | cut -f 1 | uniq
```
 
### Split vcf file by chromosome
 
```
     set -euo pipefail
     bcftools view -r ~{regions} ~{vcfFile} | bgzip -c > ~{basename}.vcf.gz
```
### Run variant effect prediction
 
```
     set -euo pipefail
 
     if [ "~{species}" = "homo_sapiens" ]; then
       human_only_command_line="--polyphen b --af --af_1kg --af_esp --af_gnomad"
     else
       human_only_command_line=""
     fi
 
     if ~{vepStats} ; then
       vepStats_command_line=""
     else 
       vepStats_command_line="--no_stats"
     fi
 
 
     vep --offline --dir ~{vepCacheDir} -i ~{vcfFile} --fasta ~{referenceFasta} --species ~{species} \
           --assembly ~{ncbiBuild} -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip ~{addParam} \
           --no_progress --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --mane \
           --canonical --protein --biotype --uniprot --tsl --variant_class --check_existing --total_length \
           --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele \
           --pick_order canonical,tsl,biotype,rank,ccds,length \
           $vepStats_command_line \
           $human_only_command_line \
           --pubmed --fork 4 --regulatory
 
```
 
### getSampleNames
 
```
     set -euo pipefail
 
     TUMR="~{tumorName}"
 
     if [ -z "~{normalName}" ]; then
         NORM="unmatched";
     else NORM="~{normalName}";
     fi
 
     echo $TUMR > names.txt
     echo $NORM >> names.txt
 
```
 
### tumorOnlyAlign task
 
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
 
### vcf2maf conversion
 
```
     set -euo pipefail
 
     TUMR=$(sed -n 1p ~{tumorNormalNames} )
     NORM=$(sed -n 2p ~{tumorNormalNames} )
 
     bgzip -c -d ~{vcfFile} > ~{basename}
 
     if ~{retainInfoProvided} ; then
         retainInfo_command_line="--retain-info MBQ,MMQ,TLOD,set"
     else
         retainInfo_command_line=""
     fi
 
     vcf2maf --ref-fasta ~{referenceFasta} --species ~{species} --ncbi-build ~{ncbiBuild} \
             --input-vcf ~{basename} --output-maf ~{basename}.maf \
             --tumor-id $TUMR --normal-id $NORM --vcf-tumor-id $TUMR --vcf-normal-id $NORM \
             --vep-path ~{vepPath} --vep-data ~{vepCacheDir} \
             --min-hom-vaf ~{minHomVaf} --buffer-size ~{bufferSize} \
             $retainInfo_command_line \
             --vep-stats ~{vepStats}
```
 
### mergeMafs
 
```
     set -euo pipefail
 
     head -n 2 ~{mafs[0]} > ~{basename}
     cat ~{sep=" " mafs} | grep -v ^# | grep -v "Hugo_Symbol" >> ~{basename}
     bgzip -c ~{basename} > ~{basename}.maf.gz

```
 
### Merging vcf files
 
```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeVcfs \
     -I ~{sep=" -I " vcfs} ~{extraArgs} \
     -O ~{basename}.vcf.gz
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

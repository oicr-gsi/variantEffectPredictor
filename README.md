# vepWorkflow

Using VEP to a vcf file and providing additional options as well

## Overview

## Dependencies

* [bedtools 2.27](https://github.com/arq5x/bedtools)
* [tabix 0.2.6](https://github.com/samtools/tabix)
* [vep 92.0](https://github.com/Ensembl/ensembl-vep)


## Usage

### Cromwell
```
java -jar cromwell.jar run vepWorkflow.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`vcfFile`|File|Input VCF file
`vcfIndex`|File|Input VCF Index File
`tumorName`|String|Name of the tumor
`normalName`|String|Name of the normal


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBed`|File?|None|Bed target file
`custom`|File?|None|Custom Appending
`toMAF`|Boolean?|None|Converting the MAF file to VEP


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBedTask.basename`|String|basename("~{vcfFile}",".vcf.gz")|
`targetBedTask.modules`|String|"bedtools/2.27 tabix/0.2.6"|Module needed to run UMI-tools extract
`targetBedTask.jobMemory`|Int|32|Memory allocated for this job (GB)
`targetBedTask.threads`|Int|4|Requested CPU threads
`targetBedTask.timeout`|Int|6|hours before task timeout
`vepAfterBed.basename`|String|basename("~{vcfFile}",".vcf.gz")|
`vepAfterBed.cacheDir`|String|"$VEP_HG19_CACHE_ROOT/.vep"|
`vepAfterBed.modules`|String|"vep/92.0 tabix/0.2.6 vep-hg19-cache/92"|
`vepAfterBed.jobMemory`|Int|32|
`vepAfterBed.threads`|Int|4|
`vepAfterBed.timeout`|Int|6|
`vcf2mafAfterBed.basename`|String|basename("~{vcfFile}",".vcf.gz")|
`vcf2mafAfterBed.modules`|String|"vcf2maf/1.6.17 tabix/0.2.6 hg19/p13 vep-hg19-cache/92 vep-hg19-exac/0.3.1"|
`vcf2mafAfterBed.referenceFasta`|String|"$HG19_ROOT/hg19_random.fa"|
`vcf2mafAfterBed.vepPath`|String|"$VEP_ROOT/bin/"|
`vcf2mafAfterBed.cacheDir`|String|"$VEP_HG19_CACHE_ROOT/.vep"|
`vcf2mafAfterBed.vcfFilter`|String|"$VEP_HG19_EXAC_ROOT/ExAC_nonTCGA.r0.3.1.somatic.sites.vep.vcf.gz"|
`vcf2mafAfterBed.jobMemory`|Int|32|
`vcf2mafAfterBed.threads`|Int|4|
`vcf2mafAfterBed.timeout`|Int|6|
`vep.basename`|String|basename("~{vcfFile}",".vcf.gz")|
`vep.cacheDir`|String|"$VEP_HG19_CACHE_ROOT/.vep"|
`vep.modules`|String|"vep/92.0 tabix/0.2.6 vep-hg19-cache/92"|Module needed to run UMI-tools extract
`vep.jobMemory`|Int|32|Memory allocated for this job (GB)
`vep.threads`|Int|4|Requested CPU threads
`vep.timeout`|Int|6|hours before task timeout
`vcf2maf.basename`|String|basename("~{vcfFile}",".vcf.gz")|base name
`vcf2maf.modules`|String|"vcf2maf/1.6.17 tabix/0.2.6 hg19/p13 vep-hg19-cache/92 vep-hg19-exac/0.3.1"|Module needed to run UMI-tools extract
`vcf2maf.referenceFasta`|String|"$HG19_ROOT/hg19_random.fa"|Reference fasta file
`vcf2maf.vepPath`|String|"$VEP_ROOT/bin/"|Path to vep script
`vcf2maf.cacheDir`|String|"$VEP_HG19_CACHE_ROOT/.vep"|Dir of cache files
`vcf2maf.vcfFilter`|String|"$VEP_HG19_EXAC_ROOT/ExAC_nonTCGA.r0.3.1.somatic.sites.vep.vcf.gz"|
`vcf2maf.jobMemory`|Int|32|Memory allocated for this job (GB)
`vcf2maf.threads`|Int|4|Requested CPU threads
`vcf2maf.timeout`|Int|6|hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`outputVcf`|File?|vcf output file
`outputMaf`|File?|optional maf output file
`outputVcfBed`|File?|Same output for vcf after bed
`outputMafBed`|File?|Same output for maf after bed


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_

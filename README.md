# variantEffectPredictor

Using VEP to a vcf file and providing additional option as well

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
`toMAF`|Boolean|If true, generate the MAF file
`onlyTumor`|Boolean|If true, run tumor only mode
`vep.vepCacheDir`|String|Directory of cache files
`vep.referenceFasta`|String|Reference fasta file
`vep.modules`|String|Required environment modules
`vcf2maf.modules`|String|Required environment modules
`vcf2maf.referenceFasta`|String|Reference fasta file
`vcf2maf.vepPath`|String|Path to vep script
`vcf2maf.vepCacheDir`|String|Directory of vep cache files
`vcf2maf.vcfFilter`|String|Filter for the vep module that is used in vcf2maf


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBed`|File?|None|Target bed file


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBedTask.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`targetBedTask.modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`targetBedTask.jobMemory`|Int|32|Memory allocated for this job (GB)
`targetBedTask.threads`|Int|4|Requested CPU threads
`targetBedTask.timeout`|Int|6|Hours before task timeout
`vep.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`vep.addParam`|String?|None|Additional vep parameters
`vep.jobMemory`|Int|32|Memory allocated for this job (GB)
`vep.threads`|Int|4|Requested CPU threads
`vep.timeout`|Int|16|Hours before task timeout
`tumorOnlyAlign.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorOnlyAlign.modules`|String|"bcftools/1.9 tabix/0.2.6 vcftools/0.1.16"|Required environment modules
`tumorOnlyAlign.jobMemory`|Int|32|Memory allocated for this job (GB)
`tumorOnlyAlign.threads`|Int|4|Requested CPU threads
`tumorOnlyAlign.timeout`|Int|6|Hours before task timeout
`getSampleNames.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`getSampleNames.modules`|String|"vcftools/0.1.16"|Required environment modules
`getSampleNames.jobMemory`|Int|32|Memory allocated for this job (GB)
`getSampleNames.threads`|Int|4|Requested CPU threads
`getSampleNames.timeout`|Int|6|Hours before task timeout
`vcf2maf.basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`vcf2maf.jobMemory`|Int|32|Memory allocated for this job (GB)
`vcf2maf.threads`|Int|4|Requested CPU threads
`vcf2maf.timeout`|Int|48|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`outputVcf`|File|Annotated vcf output file from vep
`outputTbi`|File|Index of the annotated vcf output file from vep
`outputMaf`|File?|Maf output file from vcf2maf(if toMAF is true)
`outputTargetVcf`|File?|Vcf on target for the input vcf (if targetBed is given), non annotated
`outputTargetTbi`|File?|Index of the vcf on target for the input vcf (if targetBed is given), non annotated


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

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

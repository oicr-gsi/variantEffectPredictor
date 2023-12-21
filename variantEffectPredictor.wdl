version 1.0

struct GenomeResources{
  String vep_modules
  String vcf2maf_modules
  String vepCacheDir
  String referenceFasta
  String species
  String ncbiBuild
  String vepPath
}


workflow variantEffectPredictor {
  input {
    File vcfFile
    File vcfIndex
    String? targetBed
    String tumorName
    String? normalName
    Boolean toMAF
    Boolean onlyTumor
    String reference
  }

  Map[String,GenomeResources] resources = {
    "hg19": {
       "vep_modules": "vep/105.0 tabix/0.2.6 vep-hg19-cache/105 hg19/p13",
       "vcf2maf_modules": "vcf2maf/1.6.21b tabix/0.2.6 hg19/p13 vep-hg19-cache/105",
       "vepCacheDir": "$VEP_HG19_CACHE_ROOT/.vep",
       "referenceFasta": "$HG19_ROOT/hg19_random.fa",
       "species": "homo_sapiens",
       "ncbiBuild": "GRCh37",
       "vepPath": "$VEP_ROOT/bin/"
    },
    "hg38":{
       "vep_modules": "vep/105.0 tabix/0.2.6 vep-hg38-cache/105 hg38/p12",
       "vcf2maf_modules": "vcf2maf/1.6.21b tabix/0.2.6 hg38/p12 vep-hg38-cache/105",
       "vepCacheDir": "$VEP_HG38_CACHE_ROOT/.vep",
       "referenceFasta": "$HG38_ROOT/hg38_random.fa",
       "species": "homo_sapiens",
       "ncbiBuild": "GRCh38",
       "vepPath": "$VEP_ROOT/bin/"
    },
    "mm39":{
       "vep_modules": "vep/105.0 tabix/0.2.6 vep-mm39-cache/105 mm39/p6",
       "vcf2maf_modules": "vcf2maf/1.6.21b tabix/0.2.6 mm39/p6 vep-mm39-cache/105",
       "vepCacheDir": "$VEP_MM39_CACHE_ROOT/.vep",
       "referenceFasta": "$MM39_ROOT/mm39.fa",
       "species": "mus_musculus",
       "ncbiBuild": "GRCm39",
       "vepPath": "$VEP_ROOT/bin/"
    }
  } 

  if (defined(targetBed) == true) {
    call targetBedTask {
      input: vcfFile = vcfFile,
             targetBed = targetBed
    }
  }

  if (toMAF == true) {
    call getSampleNames {
        input: tumorName = tumorName,
               normalName = normalName

    }
  }

  call chromosomeArray {
      input: vcfFile = select_first([targetBedTask.targetedVcf, vcfFile])
  }

  scatter (interval in flatten(chromosomeArray.out)) {
    call getChrCoefficient {       
      input: chromosome = interval
    }
    call subsetVcf {
      input: vcfFile = select_first([targetBedTask.targetedVcf, vcfFile]),
             vcfIndex = select_first([targetBedTask.targetedTbi, vcfIndex]),
             regions = interval,
             scaleCoefficient = getChrCoefficient.coeff
    }

    call vep {
      input: 
        vcfFile = subsetVcf.subsetVcf,
        modules = resources[reference].vep_modules,
        vepCacheDir = resources[reference].vepCacheDir,
        referenceFasta = resources[reference].referenceFasta,
        species = resources[reference].species,
        ncbiBuild = resources[reference].ncbiBuild,
        scaleCoefficient = getChrCoefficient.coeff
    }

    if (toMAF == true) {
      if (onlyTumor == true) {
        call tumorOnlyAlign {
          input: vcfFile = subsetVcf.subsetVcf,
                 tumorNormalNames = select_first([getSampleNames.tumorNormalNames]),
                 scaleCoefficient = getChrCoefficient.coeff
        }
      }
      call vcf2maf {
        input: vcfFile = select_first([tumorOnlyAlign.unmatchedOutputVcf,subsetVcf.subsetVcf]),
             tumorNormalNames = select_first([getSampleNames.tumorNormalNames]),
             vepCacheDir = resources[reference].vepCacheDir,
             modules = resources[reference].vcf2maf_modules,
             referenceFasta = resources[reference].referenceFasta,
             species = resources[reference].species,
             ncbiBuild = resources[reference].ncbiBuild,
             vepPath = resources[reference].vepPath,
             scaleCoefficient = getChrCoefficient.coeff
        }
      }
  }

  if (toMAF == true) {
    call mergeMafs {
      input: mafs = select_all(vcf2maf.mafOutput)
    }
  }

  call mergeVcfs {
    input: vcfs = vep.vepVcfOutput
  }


  parameter_meta {
    vcfFile: "Input VCF file"
    vcfIndex: "Input VCF index file"
    tumorName: "Name of the tumor sample"
    normalName: "Name of the normal sample"
    targetBed: "Target bed file"
    toMAF: "If true, generate the MAF file"
    onlyTumor: "If true, run tumor only mode"
    reference: "reference genome for input sample"
  }

  meta {
    author: "Rishi Shah, Xuemei Luo"
    email: "rshah@oicr.on.ca xuemei.luo@oicr.on.ca"
    description: "Variant Effect Predictor Workflow version 2.2"
    dependencies:
    [
      {
        name: "bedtools/2.27",
        url: "https://github.com/arq5x/bedtools"
      },
      {
        name: "tabix/0.2.6",
        url: "https://github.com/samtools/tabix"
      },
      {
        name: "vep/105.0",
        url: "https://github.com/Ensembl/ensembl-vep"
      },
      {
        name: "vcf2maf/1.6.21b",
        url: "https://github.com/mskcc/vcf2maf/commit/5ed414428046e71833f454d4b64da6c30362a89b"
      },      
      {
        name: "vcftools/0.1.16",
        url: "https://vcftools.github.io/index.html"
      }
    ]
    output_meta: {
      outputVcf: "Annotated vcf output file from vep",
      outputTbi: "Index of the annotated vcf output file from vep",
      outputMaf: "Maf output file from vcf2maf(if toMAF is true)",
      outputTargetVcf: "Vcf on target for the input vcf (if targetBed is given), non annotated",
      outputTargetTbi: "Index of the vcf on target for the input vcf (if targetBed is given), non annotated"
    }
  }

  output {
    File outputVcf = mergeVcfs.mergedVcf
    File outputTbi = mergeVcfs.mergedVcfTbi
    File? outputMaf = mergeMafs.mergedMaf
    File? outputTargetVcf = targetBedTask.targetedVcf
    File? outputTargetTbi = targetBedTask.targetedTbi

  }
}

# ================================================================
#  Scaling coefficient - use to scale RAM allocation by chromosome
# ================================================================
task getChrCoefficient {
  input {
    Int memory = 1
    Int timeout = 1
    Int largestChrom
    String chromosome
    String modules
    String refFai
  }

  parameter_meta {
    refFai: ".fai file for the reference genome, we use it to extract chromosome ids"
    timeout: "Hours before task timeout"
    chromosome: "Chromosome to check"
    memory: "Memory allocated for this job"
    modules: "Names and versions of modules to load"
    largestChrom: "Length of the largest chromosome in a genome"
  }

  command <<<
    grep -w ^~{chromosome} ~{refFai} | cut -f 2 | awk '{print int(($1/~{largestChrom} + 0.1) * 10)/10}'
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String coeff = read_string(stdout())
  }

  meta {
    output_meta: {
      coeff: "Length ratio as relative to the largest chromosome."
    }
  }
}

task targetBedTask {
  input {
    File vcfFile
    String basename = basename("~{vcfFile}", ".vcf.gz")
    File? targetBed
    String modules = "bedtools/2.27 tabix/0.2.6"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6

  }

  parameter_meta {
    vcfFile: "Vcf input files"
    targetBed: "Bed file with targets"
    basename: "Base name"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    bedtools intersect -header -u \
                       -a ~{vcfFile} \
                       -b ~{targetBed} \
                       > ~{basename}.targeted.vcf

    bgzip -c ~{basename}.targeted.vcf > ~{basename}.targeted.vcf.gz

    tabix -p vcf ~{basename}.targeted.vcf.gz
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File targetedVcf = "~{basename}.targeted.vcf.gz"
    File targetedTbi = "~{basename}.targeted.vcf.gz.tbi"
  }

  meta {
    output_meta: {
      targetedVcf: "Vcf input targeted with BED file",
      targetedTbi: "Index of the input vcf on target"
    }
  }
}


task chromosomeArray {
  input {
    File vcfFile
    Int jobMemory = 1
    Int threads = 4
    Int timeout = 1
  }

  command <<<
    zcat ~{vcfFile} | grep -v ^# | cut -f 1 | uniq
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    vcfFile: "Vcf input file"
    jobMemory: "Memory allocated to job (in GB)."
    threads: "Requested CPU threads."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }
}

task subsetVcf {
  input {
    File vcfFile
    File vcfIndex
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String regions
    String modules = "bcftools/1.9"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
    Float scaleCoefficient
  }
  command <<<
    set -euo pipefail

    bcftools view -r ~{regions} ~{vcfFile} | bgzip -c > ~{basename}.vcf.gz
  >>>

  output {
    File subsetVcf = "~{basename}.vcf.gz"
  }

  runtime {
    modules: "~{modules}"
    memory: "~{round(jobMemory * scaleCoefficient)} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    vcfFile: "Vcf input file"
    vcfIndex: "vcf index file"
    regions: "interval regions"
    basename: "Base name"
    modules: "Required environment modules"
    jobMemory: "Memory allocated to job (in GB)."
    threads: "Requested CPU threads."
    timeout: "Maximum amount of time (in hours) the task can run for."
    scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
  }
}

task vep {
  input {
    File vcfFile
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String? addParam
    String species = "homo_sapiens"
    Boolean vepStats = true
    String ncbiBuild
    String vepCacheDir
    String referenceFasta
    String modules = "vep/105.0 tabix/0.2.6 vep-hg38-cache/105 hg38/p12"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 16
    Float scaleCoefficient
  }

  parameter_meta {
    vcfFile: "Vcf input file"
    basename: "Base name"
    addParam: "Additional vep parameters"
    species: "Species name"
    ncbiBuild: "The assembly version"
    vepCacheDir: "Directory of cache files"
    referenceFasta: "Reference fasta file"
    vepStats: "If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
  }

  command <<<
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

  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{round(jobMemory * scaleCoefficient)} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File vepVcfOutput = "~{basename}.vep.vcf.gz"
  }

  meta {
    output_meta: {
      vepVcfOutput: "VEP Vcf output"
    }
  }
}

task getSampleNames {
  input {
    String tumorName
    String? normalName
    Int jobMemory = 1
    Int threads = 4
    Int timeout = 1
  }
  parameter_meta {
    tumorName: "Name of the tumor sample"
    normalName: "Name of the normal sample"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }
  command <<<
    set -euo pipefail

    TUMR="~{tumorName}"

    if [ -z "~{normalName}" ]; then
        NORM="unmatched";
    else NORM="~{normalName}";
    fi

    echo $TUMR > names.txt
    echo $NORM >> names.txt

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File tumorNormalNames = "names.txt"
  }
  meta {
    output_meta: {
      tumorNormalNames: "Names to use in the vcf2maf conversion"
    }
  }
}

task tumorOnlyAlign {
  input {
    File vcfFile
    File tumorNormalNames
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String modules = "bcftools/1.9 tabix/0.2.6"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
    Boolean updateTagValue = false
  }
  parameter_meta {
    vcfFile: "Vcf input file"
    tumorNormalNames: "Tumor and normal ID"
    basename: "Base name"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    updateTagValue: "If true, update tag values in vcf header for CC workflow"
  }

  command <<<
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
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File unmatchedOutputVcf = "~{basename}.unmatched.vcf.gz"
    File unmatchedOutputTbi = "~{basename}.unmatched.vcf.gz.tbi"
  }

  meta {
    output_meta: {
      umatchedOutputVcf: "vcf file for unmatched input",
      unmatchedOutputTbi: "index file for unmatched input"
    }
  }
}


task vcf2maf {
  input {
    File vcfFile
    String basename = basename("~{vcfFile}", ".vcf.gz")
    File tumorNormalNames
    String modules = "vcf2maf/1.6.21b tabix/0.2.6 hg38/p12 vep-hg38-cache/105"
    String species = "homo_sapiens"
    String referenceFasta
    String ncbiBuild
    String vepPath
    String vepCacheDir
    Boolean retainInfoProvided = false
    Boolean vepStats = true
    Float minHomVaf = 0.7
    Int bufferSize = 200
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 48
    Float scaleCoefficient
  }

  parameter_meta {
    vcfFile: "Vcf input file"
    species: "Species name"
    referenceFasta: "Reference fasta file"
    ncbiBuild: "The assembly version"
    vepPath: "Path to vep script"
    vepCacheDir: "Directory of vep cache files"
    retainInfoProvided: "Comma-delimited names of INFO fields to retain as extra columns in MAF"
    vepStats: "If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'"
    minHomVaf: "The minimum vaf for homozygous calls"
    bufferSize: "The buffer size"
    tumorNormalNames: "Tumor and normal ID"
    basename: "Base name"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
  }

  command <<<
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
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{round(jobMemory * scaleCoefficient)} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File mafOutput = "~{basename}.maf"
  }
  meta {
    output_meta: {
      mafOutput: "Maf output from vcf2maf"

    }
  }
}

task mergeMafs {
  input {
    Array[File] mafs
    String modules = "tabix/0.2.6"
    Int jobMemory = 24
    Int threads = 4
    Int timeout = 24
  }

  String basename = basename(mafs[0], ".maf")

  command <<<
    set -euo pipefail

    head -n 2 ~{mafs[0]} > ~{basename}
    cat ~{sep=" " mafs} | grep -v ^# | grep -v "Hugo_Symbol" >> ~{basename}
    bgzip -c ~{basename} > ~{basename}.maf.gz

  >>>

  runtime {
    modules: "~{modules}"
    memory: "~{jobMemory} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File mergedMaf = "~{basename}.maf.gz"
  }

  parameter_meta {
    mafs: "mafs from scatter to merge together."
    modules: "Required environment modules"
    jobMemory:  "Memory allocated to job (in GB)."
    threads: "Requested CPU threads."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      mergedMaf: "Merged maf"
    }
  }
}

task mergeVcfs {
  input {
    String modules = "gatk/4.1.7.0"
    Array[File] vcfs
    String? extraArgs
    Int jobMemory = 24
    Int overhead = 6
    Int threads = 4
    Int timeout = 24
  }

  String basename = basename(vcfs[0], ".vcf.gz")

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeVcfs \
    -I ~{sep=" -I " vcfs} ~{extraArgs} \
    -O ~{basename}.vcf.gz
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File mergedVcf = "~{basename}.vcf.gz"
    File mergedVcfTbi = "~{basename}.vcf.gz.tbi"
  }

  parameter_meta {
    modules: "Required environment modules."
    vcfs: "Vcf's from scatter to merge together."
    extraArgs: "Additional arguments to be passed directly to the command."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    threads: "Requested CPU threads."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      mergedVcf: "Merged vcf",
      mergedVcfTbi: "Merged vcf index"
    }
  }

}

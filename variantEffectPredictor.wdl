version 1.0 
workflow variantEffectPredictor {
  input {
    File vcfFile
    File vcfIndex
    File? targetBed
    File? custom
    Boolean? toMAF
    String tumorName
    String normalName
  }

  String customCommand = if defined(custom) == true then "true" else "false"

  if (defined(targetBed) == true) {
    call targetBedTask {
      input: vcfFile = vcfFile,
             vcfIndex = vcfIndex, 
             targetBed = targetBed
    }
    call vep as vepAfterBed {
      input: vcfFile = targetBedTask.targetedVcf,
             vcfIndex = vcfIndex,
             customCommand = customCommand,
             custom = custom
    }
    if (defined(toMAF) == true) {
      call vcf2maf as vcf2mafAfterBed {
      input: vcfFile = vepAfterBed.vepVcfOutput,
             tumorName = tumorName,
             normalName = normalName 
    }
  }
  }
  
  if (defined(targetBed) == false) {
    call vep {
      input: vcfFile = vcfFile,
             vcfIndex = vcfIndex,
             customCommand = customCommand,
	     custom = custom
    }
    if (defined(toMAF) == true) {
      call vcf2maf {
        input: vcfFile = vep.vepVcfOutput,
               tumorName = tumorName,
               normalName = normalName
      }
    }
  }
  parameter_meta {
    vcfFile: "Input VCF file"
    vcfIndex: "Input VCF Index File"
    targetBed: "Bed target file"
    custom: "Custom Appending"
    toMAF: "Converting the MAF file to VEP" 
    tumorName: "Name of the tumor"
    normalName: "Name of the normal"
  }

  meta {
    author: "Rishi Shah"
    email: "rshah@oicr.on.ca"
    description: "Using VEP to a vcf file and providing additional option as well"
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
        name: "vep/92.0",
        url: "https://github.com/Ensembl/ensembl-vep"
      }
    ]
    output_meta: {
      outputVcf: "vcf output file",
      outputMaf: "optional maf output file"
    }
  }

  output {
    File? outputVcf = vep.vepVcfOutput
    File? outputMaf = vcf2maf.mafOutput
    File? outputVcfBed = vepAfterBed.vepVcfOutput
    File? outputMafBed = vcf2mafAfterBed.mafOutput
  }
}

task targetBedTask {
  input {
    File vcfFile 
    File vcfIndex
    String basename = basename("~{vcfFile}", ".vcf.gz")
    File? targetBed
    String modules = "bedtools/2.27 tabix/0.2.6"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
  }

  parameter_meta {
    vcfFile: "vcf input files"
    vcfIndex: "Index of the vcf file"
    targetBed: "Bed file with targets"
    modules: "Module needed to run UMI-tools extract"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
  }

  command <<<
    set -euo pipefail
    bedtools intersect -header \
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
  }

  meta {
    output_meta: {
      targetedVcf: "Vcf input targeted with BED file"
    }
  }
}

task vep {
  input {
    File vcfFile 
    File vcfIndex
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String customCommand
    File? custom
    String cacheDir = "$VEP_HG19_CACHE_ROOT/.vep"
    String modules = "vep/92.0 tabix/0.2.6 vep-hg19-cache/92"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
  }

  parameter_meta {
    vcfFile: "vcf input file"
    vcfIndex: "Index of the vcf file"
    customCommand: "If the custom command is to be run"
    modules: "Module needed to run UMI-tools extract"
    custom: "Optional input for custom file"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
  }

  command <<<
    set -euo pipefail
    if [ "~{customCommand}" == "true" ]; then 
      vep --offline --cache_version 92 --dir_cache ~{cacheDir} -i ~{vcfFile} -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip --custom ~{custom} 
    fi 

    if [ "~{customCommand}" == "false" ]; then 
      vep --offline --cache_version 92 --dir_cache ~{cacheDir} -i ~{vcfFile} -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip
    fi 
  >>> 

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
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

task vcf2maf {
  input {
    File vcfFile 
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String modules = "vcf2maf/1.6.17 tabix/0.2.6 hg19/p13 vep-hg19-cache/92 vep-hg19-exac/0.3.1"
    String referenceFasta = "$HG19_ROOT/hg19_random.fa"
    String vepPath = "$VEP_ROOT/bin/"
    String cacheDir = "$VEP_HG19_CACHE_ROOT/.vep"
    String vcfFilter = "$VEP_HG19_EXAC_ROOT/ExAC_nonTCGA.r0.3.1.somatic.sites.vep.vcf.gz"
    String tumorName
    String normalName
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
  }

  parameter_meta {
    vcfFile: "vcf input file"
    referenceFasta: "Reference fasta file"
    vepPath: "Path to vep script"
    cacheDir: "Dir of cache files"
    basename: "base name"
    modules: "Module needed to run UMI-tools extract"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
  }

  command <<< 
    set -euo pipefail
    
    bgzip -c -d ~{vcfFile} > ~{basename}

    vcf2maf --ref-fasta ~{referenceFasta} --input ~{basename} --output-maf ~{basename}.maf --tumor-id ~{tumorName} --normal-id ~{normalName} --filter-vcf ~{vcfFilter} --vep-path ~{vepPath} --vep-data ~{cacheDir} 
  
    bgzip -c ~{basename}.maf > ~{basename}.maf.gz
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File mafOutput = "~{basename}.maf.gz"
  }
  meta {
    output_meta: {
      mafOutput: "maf output"
    }
  }
} 

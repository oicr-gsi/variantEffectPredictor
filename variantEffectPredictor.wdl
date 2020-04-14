version 1.0
workflow variantEffectPredictor {
  input {
    File vcfFile
    File vcfIndex
    File? targetBed
    File? custom
    Boolean toMAF
    Boolean onlyTumor
    String vepCacheModule
    String vepCacheDir
  }

  String customCommand = if defined(custom) == true then "true" else "false"

  if (defined(targetBed) == true) {
    call targetBedTask {
      input: vcfFile = vcfFile,
             vcfIndex = vcfIndex, 
             targetBed = targetBed
    }
  }
  
  call vep {
    input: vcfFile =  select_first([targetBedTask.targetedVcf, vcfFile]),
           vcfIndex = select_first([targetBedTask.targetTbi, vcfIndex]),
           customCommand = customCommand,
           custom = custom,
           vepCacheModule = vepCacheModule,
           vepCacheDir = vepCacheDir
    }
  
  if (toMAF == true) {
    if (onlyTumor == true) {
      call tumorOnlyAlign {
      input: vcfFile = vep.vepVcfOutput,
             vcfIndex = vep.vepTbiOutput
      }
    }

    call getSampleNames {
      input: vcfFile = select_first([tumorOnlyAlign.unmatchedOutputVcf, vep.vepVcfOutput]),
             vcfIndex = select_first([tumorOnlyAlign.unmatchedOutputTbi, vep.vepTbiOutput])
      }

    call vcf2maf {
      input: vcfFile = select_first([tumorOnlyAlign.unmatchedOutputVcf, vep.vepVcfOutput]),
             vcfIndex = select_first([tumorOnlyAlign.unmatchedOutputTbi, vep.vepTbiOutput]),
             vepCacheDir = vepCacheDir,
             tumorNormalNames = getSampleNames.tumorNormalNames
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
      },
      {
        name: "vcftools/0.1.16",
        url: "https://vcftools.github.io/index.html"
      },
      {
        name: "vep-hg19-cache/92"
      },
      {
        name: "vep-hg19-exac/0.3.1"
      }
    ]
    output_meta: {
      outputVcf: "vcf output file",
      outputMaf: "optional maf output file"
    }
  }

  output {
    File outputVcf = vep.vepVcfOutput
    File? outputMaf = vcf2maf.mafOutput
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
    File targetTbi = "~{basename}.targeted.vcf.gz"
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
    String vepCacheModule
    String vepCacheDir
    String modules = "vep/92.0 tabix/0.2.6 ~{vepCacheModule}"
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
      vep --offline --cache_version 92 --dir_cache ~{vepCacheDir} \
          -i ~{vcfFile} -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip --custom ~{custom} 
    fi 

    if [ "~{customCommand}" == "false" ]; then 
      vep --offline --cache_version 92 --dir_cache ~{vepCacheDir} -i ~{vcfFile} -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip
    fi 
    tabix -p vcf "~{basename}.vep.vcf.gz"
  >>> 

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File vepVcfOutput = "~{basename}.vep.vcf.gz"
    File vepTbiOutput = "~{basename}.vep.vcf.gz.tbi"

  }

  meta {
    output_meta: {
      vepVcfOutput: "VEP Vcf output",
      vepTbiOutput: "VEP Tbi output"
    }
  }
}

task tumorOnlyAlign {
  input {
    File vcfFile
    File vcfIndex
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String modules = "bcftools/1.9 vcf2maf/1.6.17 tabix/0.2.6 hg19/p13 vep-hg19-cache/92 vep-hg19-exac/0.3.1 vcftools/0.1.16"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6   
  }
  parameter_meta {
    vcfFile: "vcf input file"
    vcfIndex: "index file"
    basename: "base name"
    modules: "Module needed to run program"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
  }
  command <<<
    set -euo pipefail
    module load vcftools/0.1.16
    vcf-query -l ~{vcfFile} > sample_headers
    cat sample_headers | grep -v \"GATK\" | tr \"\\n\" \",\" > sample_names
    zcat ~{vcfFile} | sed 's/QSS\\,Number\\=A/QSS\\,Number\\=\\./' | bgzip -c > "~{basename}_input.vcf.gz"
    tabix -p vcf "~{basename}_input.vcf.gz"
    if [[ `cat sample_names | tr \",\" \"\\n\" | wc -l` == 2 ]]; then
      for item in `cat sample_names | tr \",\" \"\\n\"`; do
        if [[ $item == \"NORMAL\" || $item == *_R_* ]]; then 
          NORM=$item; else TUMR=$item;
        fi;
      done
      
    else TUMR=`cat sample_names | tr -d \",\"`; NORM=\"unmatched\"; fi
    echo -e \"$TUMR\\n$NORM\" > "~{basename}_header"
    bcftools merge ~{vcfFile} ~{vcfFile} --force-samples > "~{basename}.temp_tumor.vcf"
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
      unmatchedOutputTbi: "index file for unmatched input",
    }
  }
}

task getSampleNames {
  input {
    File vcfFile 
    File vcfIndex
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String modules = "bedtools/2.27 tabix/0.2.6 vcftools/0.1.16"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
  }
  parameter_meta {
    vcfFile: "vcf input file"
    vcfIndex: "index file"
    basename: "base name"
    modules: "Module needed to run UMI-tools extract"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
  }
  command <<<
    vcf-query -l  "~{vcfFile}" > sample_headers_all
    cat sample_headers_all | grep -v "GATK" | tr "\n" "," > sample_names_all
    if [[ `cat sample_names_all | tr "," "\n" | wc -l` == 2 ]]; then
      for item in `cat sample_names_all | tr "," "\n"`; do if [[ $item == "NORMAL" || $item == *_R_* || $item == *BC*  || $item == "unmatched" ]]; then NORM=$item; else TUMR=$item; fi; done
    else TUMR=`cat sample_names_all | tr -d ","`; NORM="unmatched"; fi
    echo $NORM > names.txt
    echo $TUMR >> names.txt
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File tumorNormalNames = "names.txt"
  }
  meta {
    output_meta: {
      tumorNormalNames: "names to use in the vcf2maf conversion"
    }
  }
}

task vcf2maf {
  input {
    File vcfFile
    File vcfIndex 
    String basename = basename("~{vcfFile}", ".vcf.gz")
    File tumorNormalNames
    String modules
    String referenceFasta
    String vepPath
    String vepCacheDir
    String vcfFilter
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
  }

  parameter_meta {
    vcfFile: "vcf input file"
    referenceFasta: "Reference fasta file"
    vepPath: "Path to vep script"
    vepCacheDir: "Dir of cache files"
    vcfFilter: "Filter for the vep module that is used in vcf2maf"
    tumorNormalNames: "Tumor and normal ID"
    basename: "base name"
    modules: "Module needed to run UMI-tools extract"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
  }

  command <<< 
    set -euo pipefail
    module unload vcftools
    NORM=$(sed -n 1p ~{tumorNormalNames} )
    TUMR=$(sed -n 2p ~{tumorNormalNames} )

    bgzip -c -d ~{vcfFile} > ~{basename}

    vcf2maf --ref-fasta ~{referenceFasta} --input ~{basename} --output-maf ~{basename}.maf \
            --tumor-id $TUMR --normal-id $NORM --vcf-tumor-id $TUMR --vcf-normal-id $NORM \
            --filter-vcf ~{vcfFilter} --vep-path ~{vepPath} --vep-data ~{vepCacheDir} 
  
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

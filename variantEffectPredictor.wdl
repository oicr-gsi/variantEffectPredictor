version 1.0
workflow variantEffectPredictor {
  input {
    File vcfFile
    File vcfIndex
    File? targetBed
    Boolean toMAF
    Boolean onlyTumor
  }

  if (defined(targetBed) == true) {
    call targetBedTask {
      input: vcfFile = vcfFile, 
             targetBed = targetBed
    }
  }

  if (toMAF == true) {
    call getSampleNames {
        input: vcfFile = vcfFile
    } 
  }

  call chromosomeArray {
      input: vcfFile = select_first([targetBedTask.targetedVcf, vcfFile])
  }

  scatter (intervals in chromosomeArray.out) {
    call subsetVcf {
      input: vcfFile = select_first([targetBedTask.targetedVcf, vcfFile]),
             vcfIndex = select_first([targetBedTask.targetedTbi, vcfIndex]),
             regions = intervals[0]
    }

    call vep {
      input: vcfFile = subsetVcf.subsetVcf
    }

    if (toMAF == true) {
      if (onlyTumor == true) {
        call tumorOnlyAlign {
          input: vcfFile = subsetVcf.subsetVcf,
                 tumorNormalNames = select_first([getSampleNames.tumorNormalNames])
        }
      }
      call vcf2maf {
        input: vcfFile = select_first([tumorOnlyAlign.unmatchedOutputVcf,subsetVcf.subsetVcf]),
             tumorNormalNames = select_first([getSampleNames.tumorNormalNames])
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
    targetBed: "Target bed file"
    toMAF: "If true, generate the MAF file"
    onlyTumor: "If true, run tumor only mode"
  }

  meta {
    author: "Rishi Shah, Xuemei Luo"
    email: "rshah@oicr.on.ca xuemei.luo@oicr.on.ca"
    description: "Variant Effect Predictor Workflow version 2.1"
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
    memory: "~{jobMemory} GB"
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
  }
}

task vep {
  input {
    File vcfFile 
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String? addParam
    String species = "homo_sapiens"
    String ncbiBuild
    String vepCacheDir
    String referenceFasta
    String modules
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 16
  }

  parameter_meta {
    vcfFile: "Vcf input file"
    basename: "Base name"
    addParam: "Additional vep parameters"
    species: "Species name"
    ncbiBuild: "The assembly version"
    vepCacheDir: "Directory of cache files"
    referenceFasta: "Reference fasta file"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
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

task getSampleNames {
  input {
    File vcfFile 
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String modules = "vcftools/0.1.16"
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 6
  }
  parameter_meta {
    vcfFile: "Vcf input file"
    basename: "Base name"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }
  command <<<
    set -euo pipefail

    vcf-query -l  "~{vcfFile}" > sample_headers_all
    cat sample_headers_all | grep -v "GATK" | tr "\n" "," > sample_names_all
    if [[ `cat sample_names_all | tr "," "\n" | wc -l` == 2 ]]; then
      for item in `cat sample_names_all | tr "," "\n"`; do if [[ $item == "NORMAL" || $item == *_R_* || $item == *_R || $item == *BC*  || $item == "unmatched" ]]; then NORM=$item; else TUMR=$item; fi; done
    else TUMR=`cat sample_names_all | tr -d ","`; NORM="unmatched"; fi

    echo $TUMR > names.txt
    echo $NORM >> names.txt

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
    String modules
    String species = "homo_sapiens"
    String referenceFasta
    String ncbiBuild
    String vepPath
    String vepCacheDir
    String vcfFilter
    Boolean retainInfoProvided = false
    Int maxfilterAC = 10
    Float minHomVaf = 0.7
    Int bufferSize = 200
    Int jobMemory = 32
    Int threads = 4
    Int timeout = 48
  }

  parameter_meta {
    vcfFile: "Vcf input file"
    species: "Species name"
    referenceFasta: "Reference fasta file"
    ncbiBuild: "The assembly version"
    vepPath: "Path to vep script"
    vepCacheDir: "Directory of vep cache files"
    vcfFilter: "Filter for the vep module that is used in vcf2maf"
    retainInfoProvided: "Comma-delimited names of INFO fields to retain as extra columns in MAF"
    maxfilterAC: "The maximum AC filter"
    minHomVaf: "The minimum vaf for homozygous calls"
    bufferSize: "The buffer size"  
    tumorNormalNames: "Tumor and normal ID"
    basename: "Base name"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<< 
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
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
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

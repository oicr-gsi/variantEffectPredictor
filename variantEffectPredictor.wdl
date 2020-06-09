version 1.0
workflow variantEffectPredictor {
  input {
    File vcfFile
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

  call vep {
      input: vcfFile =  select_first([targetBedTask.targetedVcf, vcfFile]),
  }

  if (toMAF == true) {
    if (onlyTumor == true) {
      call tumorOnlyAlign {
        input: vcfFile = select_first([targetBedTask.targetedVcf, vcfFile])    
      }
    }

    call getSampleNames {
      input: vcfFile = vcfFile
    }

    call vcf2maf {
      input: vcfFile = select_first([tumorOnlyAlign.unmatchedOutputVcf,targetBedTask.targetedVcf, vcfFile]),
             tumorNormalNames = getSampleNames.tumorNormalNames
      } 
    }

  parameter_meta {
    vcfFile: "Input VCF file"
    targetBed: "Target bed file"
    toMAF: "If true, generate the MAF file"
    onlyTumor: "If true, run tumor only mode"
  }

  meta {
    author: "Rishi Shah"
    email: "rshah@oicr.on.ca"
    description: "Variant Effect Predictor Workflow version 2.0"
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
    File outputVcf = vep.vepVcfOutput
    File outputTbi = vep.vepTbiOutput
    File? outputMaf = vcf2maf.mafOutput
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
    File targetedTbi = "~{basename}.targeted.vcf.gz.tbi"
  }

  meta {
    output_meta: {
      targetedVcf: "Vcf input targeted with BED file",
      targetedTbi: "Index of the input vcf on target"
    }
  }
}

task vep {
  input {
    File vcfFile 
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String? addParam
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
    vepCacheDir: "Directory of cache files"
    referenceFasta: "Reference fasta file"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    vep --offline --dir ~{vepCacheDir} -i ~{vcfFile} --fasta ~{referenceFasta} \
          -o ~{basename}.vep.vcf.gz --vcf --compress_output bgzip ~{addParam} \
          --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype \
          --canonical --protein --biotype --uniprot --tsl --variant_class --check_existing --total_length \
          --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele \
          --pick_order canonical,tsl,biotype,rank,ccds,length  \
          --pubmed --fork 4 --polyphen b --af --af_1kg --af_esp --af_gnomad --regulatory
    
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
    String basename = basename("~{vcfFile}", ".vcf.gz")
    String modules = "bcftools/1.9 tabix/0.2.6 vcftools/0.1.16"
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

    vcf-query -l ~{vcfFile} > sample_headers
    cat sample_headers | grep -v "GATK" | tr "\n" "," > sample_names
    zcat ~{vcfFile} | sed 's/QSS\,Number\=A/QSS\,Number\=\./' | sed 's/AS_FilterStatus\,Number\=A/AS_FilterStatus\,Number\=\./' | bgzip -c > "~{basename}_input.vcf.gz"
    tabix -p vcf "~{basename}_input.vcf.gz"
    if [[ `cat sample_names | tr "," "\n" | wc -l` == 2 ]]; then
      for item in `cat sample_names | tr "," "\n"`; do
        if [[ $item == "NORMAL" || $item == *_R_* ]]; then 
          NORM=$item; else TUMR=$item;
        fi;
      done
      
    else TUMR=`cat sample_names | tr -d ","`; NORM="unmatched"; fi
    echo -e "$TUMR\n$NORM" > "~{basename}_header"
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
      tumorNormalNames: "Names to use in the vcf2maf conversion"
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

    NORM=$(sed -n 1p ~{tumorNormalNames} )
    TUMR=$(sed -n 2p ~{tumorNormalNames} )

    bgzip -c -d ~{vcfFile} > ~{basename}

    vcf2maf --ref-fasta ~{referenceFasta} --species ~{species} --ncbi-build ~{ncbiBuild} \
            --input-vcf ~{basename} --output-maf ~{basename}.maf \
            --tumor-id $TUMR --normal-id $NORM --vcf-tumor-id $TUMR --vcf-normal-id $NORM \
            --filter-vcf ~{vcfFilter} --vep-path ~{vepPath} --vep-data ~{vepCacheDir} \
            --max-filter-ac ~{maxfilterAC} --min-hom-vaf ~{minHomVaf} --buffer-size ~{bufferSize}
  
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
      mafOutput: "Maf output from vcf2maf"

    }
  }
} 

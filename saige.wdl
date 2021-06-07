# SAIGE Workflow
# Copyright (c) 2021 Genome Research Limits
# Licensed under GPLv3, or later

# Author: Christopher Harrison <ch12@sanger.ac.uk>
# Adapted from SLURM scripts by David van Heel <d.vanheel@qmul.ac.uk>
# Incorporating Sanger requirements by Qinqin Huang <qh1@sanger.ac.uk>

version 1.0

workflow SAIGE {
  # NOTE The allosome part of the workflow has been commented out;
  # search for "No allosome yet..." in the code to reinstate.

  input {
    File plinkFile      # PLINK file for creating the GRM
    File phenoFile      # Phenotype file
    File autosomeBGENs  # File of bgen filenames, per autosomal chromosome
    File sampleFile     # File of IDs of samples in the dosage file

    # TODO No allosome yet...
    # File allosomeVCFs   # File of allosomal chromosome-VCF filename pairs, tab-delimited
  }

  # TODO Read phenotypes from file, rather than hardcoding them here
  Array[String] phenotypes = [
    # Continuous traits
    "c_random",  "c_CADSA",  "c_CADmeta", "c_BMI",     "c_HDL",
    "c_LDL",     "c_TG",

    # Binary traits
    "b_CADSA1",  "b_CADSA2", "b_HDL1",    "b_HDL2",    "b_LDL1",
    "b_LDL2",    "b_BMI1",   "b_BMI2",    "b_random1", "b_random2",
    "b_random3", "b_random4"
  ]

  # TODO Read covariants from file, rather than hardcoding them here
  Array[String] covariants = [
    "sex",  "PC1",  "PC2",  "PC3",  "PC4",  "PC5",  "PC6",
    "PC7",  "PC8",  "PC9",  "PC10", "PC11", "PC12", "PC13",
    "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20"
  ]

  # Read SPA test inputs
  Array[String] bgenFiles = read_lines(autosomeBGENs)
  # TODO No allosome yet...
  # Map[String, String] vcfFiles = read_map(allosomeVCFs)

  scatter (p in phenotypes) {
    # Step 1: Fit Null GLMM
    call FitNullGLMM {
      input:
        phenotype  = p,
        plinkFile  = plinkFile,
        phenoFile  = phenoFile,
        covariants = covariants
    }

    # Step 2: SPA Tests
    # NOTE range(22) is [0, 1, ..., 21] so we have to add 1
    scatter (chr in range(22)) {
      call SPATests_Autosome {
        input:
          phenotype     = p,
          chr           = "${chr + 1}",
          bgenFile      = bgenFiles[chr],
          sampleFile    = sampleFile,
          GMMATModel    = FitNullGLMM.GMMATModel,
          varianceRatio = FitNullGLMM.varianceRatio
      }
    }

    # TODO No allosome yet...
    # scatter (chr in ["X"]) {
    #   call SPATests_Allosome {
    #     input:
    #       phenotype     = p,
    #       chr           = chr,
    #       vcfFile       = vcfFiles[chr],
    #       vcfIndex      = "${vcfFiles[chr]}.tbi",
    #       GMMATModel    = FitNullGLMM.GMMATModel,
    #       varianceRatio = FitNullGLMM.varianceRatio
    #   }
    # }

    call Aggregate {
      input:
        phenotype    = p,
        autosomeGWAS = SPATests_Autosome.GWAS,

        # TODO No allosome yet...
        # allosomeGWAS = SPATests_Allosome.GWAS
        allosomeGWAS = []
    }
  }

  output {
    # Array of aggregated GWAS; one for each phenotype
    Array[File] GWAS = Aggregate.GWAS
  }
}

task FitNullGLMM {
  input {
    File          plinkFile
    File          phenoFile
    String        phenotype
    Array[String] covariants
  }

  command {
    # We have to do this in the shell, because there's no
    # way to extract or match against substrings in WDL
    case "${phenotype}" in
      c_*) TRAIT_TYPE="quantitative";;
      b_*) TRAIT_TYPE="binary";;
      *)   exit 1;;
    esac

    step1_fitNULLGLMM.R \
      --plinkFile=${plinkFile} \
      --phenoFile=${phenoFile} \
      --phenoCol=${phenotype} \
      --covarColList=${sep=',' covariants} \
      --sampleIDColinphenoFile=ID \
      --nThreads=8 \
      --LOCO=FALSE \
      --traitType=$TRAIT_TYPE \
      --outputPrefix=step1-${phenotype}
  }

  output {
    File GMMATModel    = "step1-${phenotype}.rda"
    File varianceRatio = "step1-${phenotype}.varianceRatio.txt"
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu:    8
    memory: "8GiB"
    zones:  "europe-west2-b" # London
  }
}

# TODO Generalise this task, so it's not split by auto- and allosome
task SPATests_Autosome {
  input {
    String phenotype
    String chr
    File   bgenFile
    File   sampleFile
    File   GMMATModel
    File   varianceRatio
  }

  command {
    step2_SPAtests.R \
      --bgenFile=${bgenFile} \
      --sampleFile=${sampleFile} \
      --minMAC=5 \
      --GMMATmodelFile=${GMMATModel} \
      --varianceRatioFile=${varianceRatio} \
      --numLinesOutput=1000 \
      --IsOutputNinCaseCtrl=TRUE \
      --IsOutputHetHomCountsinCaseCtrl=TRUE \
      --IsOutputAFinCaseCtrl=TRUE \
      --LOCO=FALSE \
      --SAIGEOutputFile=step2-${phenotype}-chr${chr}.gwas
  }

  output {
    File GWAS = "step2-${phenotype}-chr${chr}.gwas"
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu:    1
    memory: "1GiB"
    zones:  "europe-west2-b" # London
  }
}

# TODO Generalise this task, so it's not split by auto- and allosome
task SPATests_Allosome {
  input {
    String phenotype
    String chr
    File   vcfFile
    File   vcfIndex
    File   GMMATModel
    File   varianceRatio
  }

  command {
    step2_SPAtests.R \
      --vcfFile=${vcfFile} \
      --vcfFileIndex=${vcfIndex} \
      --vcfField=DS \
      --chrom=chr${chr} \
      --minMAC=5 \
      --GMMATmodelFile=${GMMATModel} \
      --varianceRatioFile=${varianceRatio} \
      --numLinesOutput=1000 \
      --IsOutputNinCaseCtrl=TRUE \
      --IsOutputHetHomCountsinCaseCtrl=TRUE \
      --IsOutputAFinCaseCtrl=TRUE \
      --LOCO=FALSE \
      --SAIGEOutputFile=step2-${phenotype}-chr${chr}.gwas
  }

  output {
    File GWAS = "step2-${phenotype}-chr${chr}.gwas"
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu:    1
    memory: "1GiB"
    zones:  "europe-west2-b" # London
  }
}

task Aggregate {
  input {
    String      phenotype
    Array[File] autosomeGWAS
    Array[File] allosomeGWAS
  }

  command <<<
    {
      # Write header
      head -1 ${autosomeGWAS[0]}

      # Write autosomal GWAS (sans header)
      tail -qn+2 ${sep=' ' autosomeGWAS}

      # TODO No allosome yet...
      # # Write allosomal GWAS (sans header, with column 3 duplicated)
      # awk 'FNR > 1 { $3 = $3 FS $3; print $0 }' ${sep=' ' allosomeGWAS}
    } \
    | gzip -c \
    > step3-${phenotype}.gwas.gz
  >>>

  output {
    File GWAS = "step3-${phenotype}.gwas.gz"
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu:    1
    memory: "1GiB"
    zones:  "europe-west2-b" # London
  }
}

# SAIGE Workflow
# Copyright (c) 2021 Genome Research Limits
# Licensed under GPLv3, or later

# Author: Christopher Harrison <ch12@sanger.ac.uk>
# Adapted from SLURM scripts by David van Heel <d.vanheel@qmul.ac.uk>
# Incorporating Sanger requirements by Qinqin Huang <qh1@sanger.ac.uk>

# NOTE This version of the workflow is for the older draft-2 WDL spec,
# using a cross product to simulate a nested scatter (as the TRE doesn't
# support subworkflows). It is thus necessarily more complicated. Please
# refer to saige.wdl for a clearer idea of what's going on.

workflow SAIGE {
  # NOTE The allosome part of the workflow has been commented out;
  # search for "No allosome yet..." in the code to reinstate.

  # TODO Read phenotypes from file, rather than hardcoding them here
  Array[String] phenotypes = [
    # Continuous Traits
    "c_random",  "c_CADSA",  "c_CADmeta", "c_BMI",     "c_HDL",
    "c_LDL",     "c_TG",

    # Binary Traits
    "b_CADSA1",  "b_CADSA2", "b_HDL1",    "b_HDL2",    "b_LDL1",
    "b_LDL2",    "b_BMI1",   "b_BMI2",    "b_random1", "b_random2",
    "b_random3", "b_random4"
  ]

  # Fit Null GLMM Inputs
  # TODO Read covariants from file, rather than hardcoding them here
  File          plinkFile  # PLINK file for creating the GRM
  File          phenoFile  # Phenotype file
  Array[String] covariants = [
    "sex",  "PC1",  "PC2",  "PC3",  "PC4",  "PC5",  "PC6",
    "PC7",  "PC8",  "PC9",  "PC10", "PC11", "PC12", "PC13",
    "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20"
  ]

  # SPA Test Inputs: Autosome
  File autosomeBGENs  # File of bgen filenames, per autosomal chromosome
  File sampleFile     # File of IDs of samples in the dosage file
  Array[String] bgenFiles = read_lines(autosomeBGENs)

  # TODO No allosome yet...
  # # SPA Test Inputs: Allosome
  # File allosomeVCFs  # File of allosomal chromosome-VCF filename pairs, tab-delimited
  # Map[String, String] vcfFiles = read_map(allosomeVCFs)

  # Step 1: Fit Null GLMM for all phenotypes
  scatter (p in phenotypes) {
    call FitNullGLMM {
      input:
        phenotype  = p,
        plinkFile  = plinkFile,
        phenoFile  = phenoFile,
        covariants = covariants
    }
  }

  # The result of a scatter is an array, which we need to index in later
  # steps. As such, we need to generate the indices manually, as we
  # cannot index by the scatter's input value.
  Array[Pair[Int, String]] enumerated_phenotypes = zip(range(length(phenotypes)), phenotypes)

  # Step 2a: SPA Tests for all phenotypes and autosomal chromosomes
  # NOTE Nested scatters aren't supported by Cromwell and subworkflows
  # aren't supported in our version of the TRE; hence the cross product
  # NOTE range(22) is [0, 1, ..., 21] so we have to add 1
  scatter (phenoXchr in cross(enumerated_phenotypes, range(22))) {
    String p     = phenoXchr.left.right
    Int    p_idx = phenoXchr.left.left
    Int    chr   = phenoXchr.right

    call SPATests_Autosome {
      input:
        phenotype     = p,
        chr           = "${chr + 1}",
        bgenFile      = bgenFiles[chr],
        sampleFile    = sampleFile,
        GMMATModel    = FitNullGLMM.GMMATModel[p_idx],
        varianceRatio = FitNullGLMM.varianceRatio[p_idx]
    }
  }

  # TODO No allosome yet...
  # # Step 2b: SPA Tests for all phenotypes and allosomal chromosomes
  # scatter (phenoXchr in cross(enumerated_phenotypes, ["X"])) {
  #   String p     = phenoXchr.left.right
  #   Int    p_idx = phenoXchr.left.left
  #   String chr   = phenoXchr.right
  #
  #   call SPATests_Allosome {
  #     input:
  #       phenotype     = p,
  #       chr           = chr,
  #       vcfFile       = vcfFiles[chr],
  #       vcfIndex      = "${vcfFiles[chr]}.tbi",
  #       GMMATModel    = FitNullGLMM.GMMATModel[p_idx],
  #       varianceRatio = FitNullGLMM.varianceRatio[p_idx]
  #   }
  # }

  # Step 3: Aggregate
  scatter (p in phenotypes) {
    call Aggregate {
      input:
        phenotype    = p,

        # These aren't used by the task directly, but are
        # needed to correctly set up the dependency graph
        autosomeGWAS = SPATests_Autosome.GWAS,

        # TODO No allosome yet...
        # allosomeGWAS = SPATests_Allosome.GWAS
        allosomeGWAS = []
    }
  }

  output {
    Array[File] GWAS = Aggregate.GWAS
  }
}

task FitNullGLMM {
  File          plinkFile
  File          phenoFile
  String        phenotype
  Array[String] covariants

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
    zones:  "europe-west2-b"  # London
  }
}

# TODO Generalise this task, so it's not split by auto- and allosome
task SPATests_Autosome {
  String phenotype
  String chr
  File   bgenFile
  File   sampleFile
  File   GMMATModel
  File   varianceRatio

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
    zones:  "europe-west2-b"  # London
  }
}

# TODO Generalise this task, so it's not split by auto- and allosome
task SPATests_Allosome {
  String phenotype
  String chr
  File   vcfFile
  File   vcfIndex
  File   GMMATModel
  File   varianceRatio

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
    zones:  "europe-west2-b"  # London
  }
}

task Aggregate {
  String phenotype

  # These aren't used by the task directly, but are
  # needed to correctly set up the dependency graph
  Array[File] autosomeGWAS
  Array[File] allosomeGWAS

  command <<<
    {
      # Autosome GWAS files
      find . -regex ".*/step2-${phenotype}-chr[1-9]+\.gwas" -type f -print0 \
      | xargs -0 awk '
        NR == 1 && FNR == 1 { print $0 }  # Print header
        FNR > 1 { print $0 }              # Print record
      '

      # Allosome GWAS files (sans header, with column 3 duplicated)
      find . -regex ".*/step2-${phenotype}-chr[XY]\.gwas" -type f -print0 \
      | xargs -0 awk 'FNR > 1 { $3 = $3 FS $3; print $0 }'
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
    zones:  "europe-west2-b"  # London
  }
}

# SAIGE Workflow
# Copyright (c) 2021 Genome Research Limits
# Licensed under GPLv3, or later

# Author: Christopher Harrison <ch12@sanger.ac.uk>
# Adapted from SLURM scripts by David van Heel <d.vanheel@qmul.ac.uk>
# Incorporating Sanger requirements by Qinqin Huang <qh1@sanger.ac.uk>

workflow SAIGE {
  File          phenotypeFile  # File of traits (one per line)
  Array[String] phenotypes = read_lines(phenotypeFile)

  # Fit Null GLMM Inputs
  String        plinkPrefix  # PLINK file prefix for creating the GRM
  File          phenoFile    # Phenotype file
  File          covarFile    # File of covariants (one per line)
  Array[String] covariants = read_lines(covarFile)

  # SPA Test Inputs
  File SPATestFOFN  # File of chromosome-SPA filename pairs, tab-delimited
  File sampleFile   # File of IDs of samples in the dosage file

  # Step 1: Fit Null GLMM for all phenotypes
  scatter (p in phenotypes) {
    call FitNullGLMM {
      input:
        phenotype   = p,
        plinkPrefix = plinkPrefix,
        phenoFile   = phenoFile,
        covariants  = covariants
    }
  }

  # The result of a scatter is an array, which we need to index in later
  # steps. As such, we need to generate the indices manually, as we
  # cannot index by the scatter's input value.
  Array[Pair[Int, String]] enumerated_phenotypes = zip(range(length(phenotypes)), phenotypes)

  # Read the SPA Test input files into a map and extract the keys
  # NOTE This relies on an undocumented feature described in
  # https://github.com/openwdl/wdl/issues/106#issuecomment-356047538
  Map[String, String] SPATestFilesByChr = read_map(SPATestFOFN)
  scatter (chr_input_pair in SPATestFilesByChr) { String chromosomes = chr_input_pair.left }

  # Step 2: SPA Tests for all phenotypes and given chromosomes
  # NOTE Nested scatters aren't supported in draft-2 WDL and
  # subworkflows aren't supported in our version of the TRE; hence the
  # cross product of arrays
  scatter (phenoXchr in cross(enumerated_phenotypes, chromosomes)) {
    String p     = phenoXchr.left.right
    Int    p_idx = phenoXchr.left.left
    String chr   = phenoXchr.right

    if (!(chr == "X" || chr == "Y")) {
      call SPATests_Autosome {
        input:
          phenotype     = p,
          chr           = chr,
          bgenFile      = SPATestFilesByChr[chr],
          sampleFile    = sampleFile,
          GMMATModel    = FitNullGLMM.GMMATModel[p_idx],
          varianceRatio = FitNullGLMM.varianceRatio[p_idx]
      }
    }

    if (chr == "X" || chr == "Y") {
      call SPATests_Allosome {
        input:
          phenotype     = p,
          chr           = chr,
          vcfFile       = SPATestFilesByChr[chr],
          vcfIndex      = "${SPATestFilesByChr[chr]}.tbi",
          GMMATModel    = FitNullGLMM.GMMATModel[p_idx],
          varianceRatio = FitNullGLMM.varianceRatio[p_idx]
      }
    }
  }

  # Step 3: Aggregate
  scatter (p in phenotypes) {
    call Aggregate {
      input:
        phenotype = p,

        # These aren't used by the task directly, but are
        # needed to correctly set up the dependency graph
        autosomeGWAS = select_all(SPATests_Autosome.GWAS),
        allosomeGWAS = select_all(SPATests_Allosome.GWAS)
    }
  }

  output {
    Array[File] GWAS = Aggregate.GWAS
  }
}

task FitNullGLMM {
  String        plinkPrefix
  File          phenoFile
  String        phenotype
  Array[String] covariants

  # These files are needed by SAIGE
  File plinkBED = "${plinkPrefix}.bed"
  File plinkBIM = "${plinkPrefix}.bim"
  File plinkFAM = "${plinkPrefix}.fam"

  command {
    # We have to do this in the shell, because there's no
    # way to extract or match against substrings in WDL
    case "${phenotype}" in
      c_*) TRAIT_TYPE="quantitative";;
      b_*) TRAIT_TYPE="binary";;
      *)   exit 1;;
    esac

    # Strip the extension off the TRE-transformed
    # PLINK files to get the actual prefix
    PLINK_FILE=$(printf "%s" ${plinkBED} | sed "s/\.bed$//")

    step1_fitNULLGLMM.R \
      --plinkFile=$PLINK_FILE \
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

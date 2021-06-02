# SAIGE Binary Traits Workflow
# Copyright (c) 2021 Genome Research Limits
# Licensed under GPLv3, or later

# Author: Christopher Harrison <ch12@sanger.ac.uk>
# Adapted from SLURM scripts by David van Heel <d.vanheel@qmul.ac.uk>

workflow SAIGE_BinaryTraits {
  # TODO Read phenotypes from file, rather than hardcoding them here
  Array[String] phenotypes = [
    "ANA_A2_V2",
    "ANA_B1_V2",
    "ANA_B2_V2",
    "ANA_C2_v2"
  ]

  # Fit Null GLMM Inputs
  File plinkFile  # PLINK file for creating the GRM
  File phenoFile  # Phenotype file

  # SPA Test Inputs: Autosome
  File autosomeBGENs  # File of bgen filenames, per autosomal chromosome
  File sampleFile     # File of IDs of samples in the dosage file
  Array[String] bgenFiles = read_lines(autosomeBGENs)

  # SPA Test Inputs: Allosome
  File allosomeVCFs  # File of allosomal chromosome-VCF filename pairs, tab-delimited
  Map[String, String] vcfFiles = read_map(allosomeVCFs)

  scatter (p in phenotypes) {
    # Step 1: Fit Null GLMM
    call FitNullGLMM {
      input:
        phenotype = p,
        plinkFile = plinkFile,
        phenoFile = phenoFile
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

    scatter (chr in ["X"]) {
      call SPATests_Allosome {
        input:
          phenotype     = p,
          chr           = chr,
          vcfFile       = vcfFiles[chr],
          vcfIndex      = "${vcfFiles[chr]}.tbi",
          GMMATModel    = FitNullGLMM.GMMATModel,
          varianceRatio = FitNullGLMM.varianceRatio
      }
    }
  }

  output {
    # TODO
  }
}

task FitNullGLMM {
  File   plinkFile
  File   phenoFile
  String phenotype

  # TODO Should the list of covariants be an input?
  command {
    step1_fitNULLGLMM.R \
      --plinkFile=${plinkFile} \
      --phenoFile=${phenoFile} \
      --phenoCol=${phenotype} \
      --covarColList=EthnicTightPCAPakBang_Dec2020,age_at_recruit,age_at_recruit_sq,S1QST_Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
      --sampleIDColinphenoFile=SAMPLE_ID \
      --nThreads=8 \
      --LOCO=FALSE \
      --traitType=binary \
      --outputPrefix=step1-${phenotype}
  }

  output {
    File GMMATModel    = "step1-${phenotype}.rda"
    File varianceRatio = "step1-${phenotype}.varianceRatio.txt"
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu: 8
    memory: "8 GB"

    # London
    zones: "europe-west2-b"
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
    cpu: 1
    memory: "1 GB"

    # London
    zones: "europe-west2-b"
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
    cpu: 1
    memory: "1 GB"

    # London
    zones: "europe-west2-b"
  }
}

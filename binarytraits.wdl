# SAIGE Binary Traits Workflow
# Copyright (c) 2021 Genome Research Limits
# Licensed under GPLv3, or later

# Author: Christopher Harrison <ch12@sanger.ac.uk>
# Adapted from SLURM scripts by David van Heel <d.vanheel@qmul.ac.uk>

workflow SAIGE_BinaryTraits {
  File plinkFile
  File phenoFile

  Array[String] phenotypes = [
    "ANA_A2_V2",
    "ANA_B1_V2",
    "ANA_B2_V2",
    "ANA_C2_v2"
  ]

  scatter (pheno in phenotypes) {
    call FitNullGLMM {
      input:
        plinkFile = plinkFile,
        phenoFile = phenoFile,
        phenoCol  = pheno
    }
  }

  output {
    # TODO
  }
}

task FitNullGLMM {
  File   plinkFile
  File   phenoFile
  String phenoCol

  command {
    step1_fitNULLGLMM.R \
      --plinkFile=${plinkFile} \
      --phenoFile=${phenoFile} \
      --phenoCol=${phenoCol} \
      --covarColList=EthnicTightPCAPakBang_Dec2020,age_at_recruit,age_at_recruit_sq,S1QST_Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
      --sampleIDColinphenoFile=SAMPLE_ID \
      --nThreads=8 \
      --LOCO=FALSE \
      --traitType=binary \
      --outputPrefix=step1-${phenoCol}
  }

  output {
    File rda            = "step1-${phenoCol}.rda"
    File varianceRation = "step1-${phenoCol}.varianceRatio.txt"
    # TODO Other outputs?
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu: 8
    memory: "8 GB"

    # London
    zones: "europe-west2-b"
  }
}

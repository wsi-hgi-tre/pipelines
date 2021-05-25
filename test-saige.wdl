workflow test_saige_pipeline {
  call test_saige

  output {
    test_saige.out
  }
}

task test_saige {
  command {
    step1_fitNULLGLMM.R --help
  }

  output {
    String out = read_string(stdout())
  }

  runtime {
    docker: "eu.gcr.io/fg-qmul-testing-master/saige:0.44.5"
    cpu: 1
    memory: "8 GB"

    # London
    zones: "europe-west2-b"
  }
}

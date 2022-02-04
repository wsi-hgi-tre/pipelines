workflow test_container_prspipe {                                                                                                     
  String input_string                                                                                                         
  call example_task {input: input_string = input_string}                                                                      
  output { Array[String] message = example_task.message }                                                                     
}                                                                                                                             
                                                                                                                              
task example_task {                                                                                                           
  String input_string                                                                                                         
                                                                                                                              
  command {                                                                                                                   
    echo "testing container. Input string: ${input_string}"                                                                   
    /usr/local/bin/plink2 --version                                                                                           
    /usr/local/bin/plink2_x86_64 --version                                                                                    
    /usr/local/bin/plink --version                                                                                            
    /usr/local/bin/qctool --version                                                                                           
  }                                                                                                                           
                                                                                                                              
  output {                                                                                                                    
  Array[String] message = read_lines(stdout())                                                                                
  }                                                                                                                           
                                                                                                                              
  runtime {                                                                                                                   
  docker: "eu.gcr.io/fg-qmul-containers/wtsihgi/tre_prspipe:d17ea65"                                                          
  cpu: 1                                                                                                                      
  memory: "1GiB"                                                                                                              
  zones: "europe-west2-b" # london                                                                                            
  }                                                                                                                           
                                                                                                                              
}                                                                                                                             
     

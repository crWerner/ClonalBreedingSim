

if (dominance_degree > 0) {  # see "0_GlobalParameters.R"
  
  seed_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 15))
  colnames(seed_output) =c("cycle", "scenario", "stage", "gain", "meanA", "meanD", "varG", 
                          "varA", "varD", "AvsD" ,"h2", "H2", "rEBV", "rEGV", "Fcoeff")
  seed_output$stage = "SEED"
  
  st1_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 15))
  colnames(st1_output) =c("cycle", "scenario", "stage", "gain", "meanA", "meanD", "varG", 
                         "varA", "varD", "AvsD" , "h2", "H2", "rEBV", "rEGV", "Fcoeff")
  st1_output$stage = "ST1"
  
  st5_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 15))
  colnames(st5_output) =c("cycle", "scenario", "stage", "gain", "meanA", "meanD", "varG", 
                         "varA", "varD", "AvsD" , "h2", "H2", "rEBV", "rEGV", "Fcoeff")
  st5_output$stage = "ST5"
  
} else {
  
  seed_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 9))
  colnames(seed_output) <- c("cycle", "scenario", "stage", "gain", "meanA", 
                             "varG", "H2", "rEBV", "Fcoeff")
  seed_output$stage = "SEED"
  
  st1_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 9))
  colnames(st1_output) <- c("cycle", "scenario", "stage", "gain", "meanA",  
                            "varG", "H2", "rEBV", "Fcoeff")
  st1_output$stage = "ST1"
  
  st5_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 9))
  colnames(st5_output) <- c("cycle", "scenario", "stage", "gain", "meanA", 
                            "varG", "H2", "rEBV", "Fcoeff")
  st5_output$stage = "ST5"
  
}









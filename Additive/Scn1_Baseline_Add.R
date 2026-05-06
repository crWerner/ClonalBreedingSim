
rm(list = ls())

library(AlphaSimR)
load("Burnin_RTB_additive.RData")

######################### BASELINE SCENARIO ###################################

cat(paste0("\n","Phenotypic Selection","\n"))

for(c in 1:cycles){
  
  b <- brnIn + c 
  cat(paste0("cycle ", b, "\n"))
  
  cat(" -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_Phe_Add.R")

  cat(" -> Update Parents\n")   
  source("2-1_UpdateParents_Phe_Add.R")

  cat(" -> Advance Year\n")
  source("2-2_AdvanceYear_Phe_Add.R")
  
}


### CALCULATE THE POPULATION PARAMETERS FOR THE FINAL POP

b <- b + 1

cat(" -> Calculate Population Parameters\n")
popPar_seed <- genParam(seedlings)
popPar_St1 <- genParam(ST1)
popPar_St5 <- genParam(ST5)

source("PopParam_Phe_Add.R")

seed_output$scenario <- st1_output$scenario <- st5_output$scenario <- "Pheno"
output <- rbind(seed_output, st1_output, st5_output)

saveRDS(output,paste0("Scn1_Baseline_add.rds"))


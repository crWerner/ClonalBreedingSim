
rm(list=ls())

library(AlphaSimR)
load("BurnIn_RTB_additive.RData")


######################### GS in Stage1 SCENARIO ################################

n_famLines <- n_famLinesGS  # REDUCE NUMBER OF INDIVIDUALS PER CROSS

cat(paste0("\n", "GS in Stage 1", "\n"))

for(c in 1:cycles){ 

  b <- brnIn + c 
  cat(paste0("cycle ", b, "\n"))
  
  cat(" -> Calculate EBVs and EGVs\n")
  
  train_pop <- train_pop[(nInd(train_pop) - 1501):nInd(train_pop)] # for testing script
  
  gs_model <- RRBLUP(pop = train_pop)     
  
  if (c == 1) parents <- setEBV(parents, gs_model)
  ST1 <- setEBV(ST1, gs_model)
  
  cat("  -> Update Parents\n")
  source("3-1_UpdateParents_GS-St1_Add.R")

  
  cat(" -> Calculate Population Parameters\n")   
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_GS-St1_Add.R")      # calculation of EGV for ST1 within script

  
  fx_eff <- b
  cat("  -> Advance Year\n")
  source("3-2_AdvanceYear_GS-St1_Add.R")
  
  # New seedlings - after Advance Year to use same script 3_2 for normal ST1-GS and ST-GS-OCS
  seedlings <- randCross(parents, nCrosses = n_crosses, nProgeny = n_famLines)
  seedlings <- setPheno(seedlings, varE = varE_seed, reps = 1, fixEff = fx_eff)
  
  
  cat("  -> Train Prediction Model\n")
  
  if (c == 1) {
    train_pop <- mergePops(list(train_pop, ST1, ST2, ST3, ST4))
    
  } else {
    train_pop <- mergePops(list(train_pop, ST1, ST2, ST3, ST4, ST5))
    
  }
}


### CALCULATE THE POPULATION PARAMETERS FOR THE FINAL POP

b <- b + 1

train_pop <- train_pop[(nInd(train_pop) - 1501):nInd(train_pop)] # for testing script
gs_model <- RRBLUP(train_pop)     

ST1 <- setEBV(ST1, gs_model)

cat(" -> Calculate Population Parameters\n")
popPar_seed <- genParam(seedlings)
popPar_St1 <- genParam(ST1)
popPar_St5 <- genParam(ST5)

source("PopParam_GS-St1_Add.R")   # calculation of EGV for ST1 within script

seed_output$scenario <- st1_output$scenario <- st5_output$scenario <- "ST1GS"

output <- rbind(seed_output, st1_output, st5_output)

saveRDS(output, "Scn2_GS-St1_add.rds")


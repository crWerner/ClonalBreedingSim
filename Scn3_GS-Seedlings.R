
rm(list = ls())

library(AlphaSimR)
load("BurnIn_RTB.RData")


######################### GS in seedlings SCENARIO #############################

n_famLines <- n_famLinesGSseed
nCrosses <-  n_crosses_GSseed


cat(paste0("\n", "GS in seedlings", "\n"))

for (c in 1:cycles) { 
  
  b <- brnIn + c
  cat(paste0("cycle ", b, "\n"))


  cat(" -> Calculate EBVs and EGVs\n")
  
  train_pop <- train_pop[(nInd(train_pop) - 1501):nInd(train_pop)] # for testing script
  
  if (dominance_degree > 0) {
    
    gs_model <- RRBLUP_D(pop = train_pop, maxIter = 80)     
    
    # "bv" = genomic estimated breeding values (GEBV)
    seedlings <- setEBV(seedlings, gs_model, value = "bv")
    ST1 <- setEBV(ST1, gs_model, value = "bv")
    
  } else {
    
    gs_model <- RRBLUP(pop = train_pop)     
    
    seedlings <- setEBV(seedlings, gs_model)
    ST1 <- setEBV(ST1, gs_model)
    
  }
  
  cat("  -> Update Parents\n")
  source("4-1_UpdateParents_GS-seedlings.R")
  
  
  cat("  -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_GS_Seedlings.R")
  
  fx_eff <- b
  cat(" -> Advance Year\n")
  source("4-2_AdvanceYear_GS-seedlings.R")
  
  # New seedlings - after Advance Year to use same script 4_2 for normal ST1-GS and ST-GS-OCS
  seedlings <- randCross(parents, nCrosses = n_crosses, nProgeny = n_famLines)
  seedlings <- setPheno(seedlings, varE = varE_seed, reps = 1, fixEff = fx_eff)
  
  
  cat(" -> Train Prediction Model\n")
  
  if (c == 1) {
    train_pop <- mergePops(list(train_pop, ST1, ST2, ST3, ST4))
    
  } else {
    train_pop <- mergePops(list(train_pop, ST1, ST2, ST3, ST4, ST5))
    
  }
}

### CALCULATE THE POPULATION PARAMETERS FOR THE FINAL POP

b <- b + 1

train_pop <- train_pop[(nInd(train_pop) - 1501):nInd(train_pop)] # for testing script

if (dominance_degree > 0) {
  
  gs_model <- RRBLUP_D(pop = train_pop, maxIter = 80)     
  # "bv" = genomic estimated breeding values (GEBV)
  seedlings <- setEBV(seedlings, gs_model, value = "bv")
  ST1 <- setEBV(ST1, gs_model, value = "bv")
  
} else {
  
  gs_model <- RRBLUP(pop = train_pop)     
  seedlings <- setEBV(seedlings, gs_model)
  ST1 <- setEBV(ST1, gs_model)
  
}

cat(" -> Calculate Population Parameters\n")
popPar_seed <- genParam(seedlings)
popPar_St1 <- genParam(ST1)
popPar_St5 <- genParam(ST5)

source("PopParam_GS_Seedlings.R")

seed_output$scenario <- st1_output$scenario <- st5_output$scenario <- "Seed_GS"

output <- rbind(seed_output, st1_output, st5_output)

saveRDS(output, "Scn3_GS-Seedlings.rds")

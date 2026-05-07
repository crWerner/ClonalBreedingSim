
rm(list = ls())

library(AlphaSimR)
load("BurnIn_RTB.RData")


######################### 2-part GS in Seedlings SCENARIO ######################

n_famLines <- n_famLinesGSseed2
n_crosses <- n_crosses_GSseed

f1_acc <- matrix(NA, nrow = cycles + 1, ncol = 3) # 3 cycles per year

cat(paste0("\n", "2-part GS in seedlings", "\n"))

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
 
  cat(" -> Update Parents\n")
  f1_acc[c, 1] <- cor(seedlings@ebv, seedlings@gv) # seedling gebv accuracy
  source("4-1_UpdateParents_GS-seedlings.R")
  
  
  cat(" -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_GS_Seedlings.R")
  
  fx_eff <- b
  cat(" -> Advance Year\n")
  source("4-2_AdvanceYear_GS-seedlings.R")
  
  
  ### New seedlings
  
  # 1st crossing cycle
  seedlings <- randCross(parents, nCrosses = n_crosses, nProgeny = n_famLines)
  
  # 2nd crossing cycle
  if (dominance_degree > 0) {
    seedlings <- setEBV(seedlings, gs_model, value = "bv")
  } else {
    seedlings <- setEBV(seedlings, gs_model)
  }
  
  f1_acc[c, 2] <- cor(seedlings@ebv, seedlings@gv)
  parents <- selectInd(seedlings, nInd = n_parents, use = "ebv")   
  seedlings <- randCross(parents, nCrosses = n_crosses, nProgeny = n_famLines)
  
  
  # 3rd crossing cycle
  if (dominance_degree > 0) {
    seedlings <- setEBV(seedlings, gs_model, value = "bv")
  } else {
    seedlings <- setEBV(seedlings, gs_model)
  }
  
  f1_acc[c, 3] <- cor(seedlings@ebv, seedlings@gv)
  parents <- selectInd(seedlings, nInd = n_parents, use = "ebv")   
  seedlings <- randCross(parents, nCrosses = n_crosses, nProgeny = n_famLines)
  
  
  cat(" -> Train Prediction Model\n")
  
  if (c == 1) {
    train_pop <- mergePops(list(train_pop, ST1, ST2, ST3, ST4))
    
  } else {
    train_pop <- mergePops(list(train_pop, ST1, ST2, ST3, ST4, ST5))
    
  }
}

### CALCULATE THE POPULATION PARAMETERS FOR THE FINAL POP

b <- b + 1

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

f1_acc[c + 1, 1] <- cor(seedlings@ebv, seedlings@gv) 

cat("  -> Calculate Population Parameters\n")
popPar_seed <- genParam(seedlings)
popPar_St1 <- genParam(ST1)
popPar_St5 <- genParam(ST5)

source("PopParam_GS_Seedlings.R")

seed_output$scenario <- st1_output$scenario <- st5_output$scenario <- "Seed_GS-3cyc"

output <- list(rbind(seed_output, st1_output, st5_output))
output[[2]] <- f1_acc

saveRDS(output, "Scn4_GS-seedlings-2_add.rds")

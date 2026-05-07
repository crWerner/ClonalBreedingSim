
rm(list = ls())

library(AlphaSimR)

load("BurnIn_RTB.RData")

### GPCP algorithm
if (dominance_degree > 0) {
  Rcpp::sourceCpp("GPCP_Dom.cpp")
} else {
  Rcpp::sourceCpp("GPCP_Add.cpp")
}


######################### GPCP in Seedlings SCENARIO ###########################

n_famLines <- n_famLinesGSseed
nCrosses <-  n_crosses_GSseed

ocs_degree <- 0 # i.e., No OCS
target_angle <- ocs_degree * pi / 180 # OCS angle in rad

n_parents_GPCP <- matrix(NA, nrow = cycles + 1, ncol = 1)


cat(paste0("\n", "GPCP in seedlings", "\n"))

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
 
  # All seedlings are potential new parents. However, since the GA would be 
  # overwhelmed with that many seedling, a within-family pre-selection of the 
  # N best individuals (highest EBV) is made first. This is for illustration only.

  # parents <- selectWithinFam(seedlings, fam_maxSeed_GPCP, use = "ebv") 

  cat("  -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_GS_Seedlings.R")
  
  fx_eff <- b
  cat(" -> Advance Year\n")
  source("4-2_AdvanceYear_GS-seedlings.R")
  
  cat(" -> Create Crossing Plan\n")
  # Identify the next generation's parents and their cross combinations

  if (dominance_degree > 0) {
    
    plan <- selectCrossPlan(
      cycleNumber = b,
      nCross = n_crosses,
      M = pullSnpGeno(parents),
      a = gs_model@gv[[1]]@addEff,
      d = gs_model@gv[[1]]@domEff, 
      targetAngle = target_angle)  
    
  } else {
    
    plan <- selectCrossPlan(
      cycleNumber = b,
      nCross = n_crosses,
      M = pullSnpGeno(parents),
      a = gs_model@gv[[1]]@addEff,
      targetAngle = target_angle)
    
  }
  
  n_parents_GPCP[c] <- length(unique(as.vector(plan$crossPlan)))
  
  cat(" -> Create Seedlings\n")
  seedlings <- makeCross(parents, plan$crossPlan, nProgeny = n_famLines) 
  seedlings <- setPheno(seedlings, varE = varE_seed, reps = 1, fixEff=fx_eff)
  
  
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


fx_eff <- b
source("4-2_AdvanceYear_GS-seedlings.R")


if (dominance_degree > 0) {
  plan <- selectCrossPlan(
    cycleNumber = b,
    nCross = n_crosses,
    M = pullSnpGeno(parents),
    a = gs_model@gv[[1]]@addEff,
    d = gs_model@gv[[1]]@domEff, 
    targetAngle = target_angle)  
  
} else {
  plan <- selectCrossPlan(
    cycleNumber = b,
    nCross = n_crosses,
    M = pullSnpGeno(parents),
    a = gs_model@gv[[1]]@addEff,
    targetAngle = target_angle)
}

n_parents_GPCP[c] <- length(unique(as.vector(plan$crossPlan)))

seed_output$scenario <- st1_output$scenario <- st5_output$scenario <- "Seed_GPCP"

output <- list(rbind(seed_output, st1_output, st5_output))
output[[2]] <- n_parents_GPCP

saveRDS(output, "Scn6_GPGCP-Seedlings.rds")

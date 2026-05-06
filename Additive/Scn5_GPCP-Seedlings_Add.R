
rm(list = ls())

library(AlphaSimR)
Rcpp::sourceCpp("GPCP_Add.cpp")
load("BurnIn_RTB_additive.RData")

# ????? source("makeCrossMultProgeny.R") 

######################### GS in Seedlings SCENARIO #############################

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
  gs_model <- RRBLUP(pop = train_pop)  
  
  seedlings <- setEBV(seedlings, gs_model)
  ST1 <- setEBV(ST1, gs_model)


 
  cat("  -> Update Parents\n")
 
   # All seedlings are poptential new parents. However, since the GA would be 
  # overwhelmed with that many seedling, a within-family pre-selection of the 
  # N best individuals (highest EBV) is made first.

  parents <- selectWithinFam(seedlings, fam_maxSeed_GPCP, use = "ebv") 

  
  cat("  -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_GS_Seedlings_Add.R")
  
  fx_eff <- b
  cat(" -> Advance Year\n")
  
  source("4-2_AdvanceYear_GS-seedlings_Add.R")
  
  cat(" -> Create Crossing Plan\n")
  # Identify the next generation's parents and their cross combinations

  plan <- selectCrossPlan(
    cycleNumber = b,
    nCross = n_crosses,
    M = pullSnpGeno(parents),
    a = gs_model@gv[[1]]@addEff,
    targetAngle = target_angle)
  
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
gs_model <- RRBLUP(train_pop)   

seedlings <- setEBV(seedlings, gs_model)
ST1 <- setEBV(ST1, gs_model)

cat(" -> Calculate Population Parameters\n")
popPar_seed <- genParam(seedlings)
popPar_St1 <- genParam(ST1)
popPar_St5 <- genParam(ST5)
source("PopParam_GS_Seedlings_Add.R")


fx_eff <- b
source("4-2_AdvanceYear_GS-seedlings_Add.R")

plan <- selectCrossPlan(
  cycleNumber = b,
  nCross = n_crosses,
  M = pullSnpGeno(parents),
  a = gs_model@gv[[1]]@addEff,
  targetAngle = target_angle)

n_parents_GPCP[c] <- length(unique(as.vector(plan$crossPlan)))


seed_output$scenario <- st1_output$scenario <- st5_output$scenario <- "Seed_GPCP"

output <- list(rbind(seed_output, st1_output, st5_output))
output[[2]] <- n_parents_GPCP

saveRDS(output, "Scn5_GPGCP-Seedlings_add.rds")

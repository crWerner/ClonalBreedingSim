
### 1. LOAD PACKAGES AND FUNCTIONS

library(AlphaSimR)
source("0_GlobalParameters_Add.R")


### GENERATE FOUNDER POPULATION

FOUNDERPOP <- runMacs(nInd = n_founders, 
                      nChr = n_chr, 
                      segSites = n_qtl + n_snp,
                      inbred = FALSE, 
                      species = species)


SP <- SimParam$new(FOUNDERPOP)
SP$addSnpChip(nSnpPerChr = n_snp)  # add SNP array

SP$addTraitAG(nQtlPerChr = n_qtl,
              mean = init_meanG,
              var = init_varG)

founderpop <- newPop(FOUNDERPOP) # create simulation object

varE_seed <- var_e(var_g = varG(founderpop), H2 = H2_seed) # varE for target H2 in seedlings
varE_St1 <- var_e(var_g = varG(founderpop), H2 = H2_St1)   # varE for target H2 in Stage 1 trials

parents <- setPheno(founderpop, varE = varE_St1, reps = rep_S3)   

source("1_FillPipeline_Add.R")


### CREATE OUTPUT MATRICES

seed_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 10))
colnames(seed_output) <- c("cycle", "scenario", "stage", "gain", "meanA","meanG", 
                         "varG", "h2", "rEBV", "Fcoeff")
seed_output$stage = "SEED"

st1_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 10))
colnames(st1_output) <- c("cycle", "scenario", "stage", "gain", "meanA", "meanG", 
                       "varG", "h2", "rEBV", "Fcoeff")
st1_output$stage = "ST1"

st5_output <- data.frame(matrix(NA, nrow = brnIn + cycles + 1, ncol = 10))
colnames(st5_output) <- c("cycle", "scenario", "stage", "gain", "meanA", "meanG", 
                       "varG", "h2", "rEBV", "Fcoeff")
st5_output$stage = "ST5"



### MEASURE INITIAL PARAMETERS OF THE CREATED BREEDING PROGRAM STAGES
### SEEDLINGS, ST1, ST2 AND ST5

# Calculate initial genomic heterozygosity of the population

seed_het <- pullQtlGeno(seedlings)
seed_het_0 <- mean(rowMeans(1 - abs(seed_het - 1)))

st1_het <- pullQtlGeno(ST1)
st1_het_0 <- mean(rowMeans(1 - abs(st1_het - 1)))

st5_het <- pullQtlGeno(ST5)
st5_het_0 <- mean(rowMeans(1 - abs(st5_het - 1)))



#### BURN IN START ########################

cat(paste0("\n","Burn In","\n"))

for(b in 1:brnIn){
  
  cat(paste0("cycle ",b,"\n"))

  cat("  -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_Phe_Add.R")
  
  cat("  -> Update Parents\n")
  source("2-1_UpdateParents_Phe_Add.R")
 

  fx_eff <- b                     # fixed effect for year
  cat("  -> Advance Year\n")
  source("2-2_AdvanceYear_Phe_Add.R")
  
  
  # Training population populated over last 3 burn-in years
  
  if(b == (brnIn - 3)){
    cat("  -> Train Prediction Model\n")
    train_pop <- ST1
  }
  if(b > (brnIn - 2)){
    train_pop <- mergePops(list(train_pop, ST1, ST2))
  }
  if(b > (brnIn - 1)){
    train_pop <- mergePops(list(train_pop, ST3))
  }
  
}

save.image("BurnIn_RTB_additive.RData")



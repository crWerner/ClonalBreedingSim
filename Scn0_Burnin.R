
### 1. LOAD PACKAGES AND FUNCTIONS

library(AlphaSimR)
source("0_GlobalParameters.R")


### GENERATE FOUNDER POPULATION

FOUNDERPOP <- runMacs(nInd = n_founders, 
                      nChr = n_chr, 
                      segSites = n_qtl + n_snp,
                      inbred = FALSE, 
                      species = species)


SP <- SimParam$new(FOUNDERPOP)

if (dominance_degree > 0) {  # see "0_GlobalParameters.R"
  
  SP$addTraitADG(nQtlPerChr = n_qtl,
                 mean = init_meanG,
                 var = init_varG,
                 meanDD = dominance_degree, 
                 varDD = dominance_degree_var, 
                 useVarA = TRUE)
  
} else {
  
  SP$addTraitAG(nQtlPerChr = n_qtl,
                mean = init_meanG,
                var = init_varG)
  
}

SP$addSnpChip(nSnpPerChr = n_snp)  # add SNP array

founderpop <- newPop(FOUNDERPOP) # create simulation object

# Calculate error variance for target H2 in seedlings and Stage 1 trials.
varE_seed <- var_e(var_g = varG(founderpop), H2 = H2_seed) 
varE_St1 <- var_e(var_g = varG(founderpop), H2 = H2_St1)   

parents <- setPheno(founderpop, varE = varE_St1, reps = rep_S3)   

source("1_FillPipeline.R")

source("output_matrix.R") # Matrices to store output data.


### Calculate initial genomic heterozygosity of the population

seed_het <- pullQtlGeno(seedlings)
seed_het_0 <- mean(rowMeans(1 - abs(seed_het - 1)))

st1_het <- pullQtlGeno(ST1)
st1_het_0 <- mean(rowMeans(1 - abs(st1_het - 1)))

st5_het <- pullQtlGeno(ST5)
st5_het_0 <- mean(rowMeans(1 - abs(st5_het - 1)))



#### BURN-IN START ########################

cat(paste0("\n","Burn In","\n"))

for(b in 1:brnIn){
  
  cat(paste0("cycle ", b, "\n"))

  cat(" -> Calculate Population Parameters\n")
  popPar_seed <- genParam(seedlings)
  popPar_St1 <- genParam(ST1)
  popPar_St5 <- genParam(ST5)
  
  source("PopParam_Phe.R")
  
  cat("  -> Update Parents\n")
  source("2-1_UpdateParents_Phe.R")
 

  fx_eff <- b                     # fixed effect for year
  cat("  -> Advance Year\n")
  source("2-2_AdvanceYear_Phe.R")
  
  
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

save.image("BurnIn_RTB.RData")

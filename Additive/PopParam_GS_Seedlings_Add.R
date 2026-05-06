### CALCULATE POPULATION PARAMETERS FOR SEEDLINGS, ST1, ST2 AND ST5

# Genetic mean, genetic variance and initial genomic heterozygosity (fraction of heterozygous QTL amongst all loci)

seed_output$cycle[b] <- b
seed_output$gain[b] <- meanG(seedlings)     
seed_output$meanA[b] <- seed_output$meanG[b] <- mean(popPar_seed$gv_a)
seed_output$varG[b] <- varG(seedlings)  
seed_output$h2[b] <- varG(seedlings) / varP(seedlings) 

seed_het <- pullQtlGeno(seedlings)
seed_het <- mean(rowMeans(1 - abs(seed_het - 1)))

seed_output$Fcoeff[b] <- 1 - (seed_het / seed_het_0)
seed_output$rEBV[b] <- cor(seedlings@ebv, seedlings@gv)


st1_output$cycle[b] <- b
st1_output$gain[b] <- meanG(ST1)
st1_output$meanA[b] <- st1_output$meanG[b] <- mean(popPar_St1$gv_a)
st1_output$varG[b] <- varG(ST1)
st1_output$h2[b] <- varG(ST1) / varP(ST1) 

st1_het <- pullQtlGeno(ST1)
st1_het <- mean(rowMeans(1 - abs(st1_het - 1)))

st1_output$Fcoeff[b] <- 1 - (st1_het / st1_het_0)
st1_output$rEBV[b] <- cor(ST1@ebv, ST1@gv)


st5_output$cycle[b] <- b
st5_output$gain[b] <- meanG(ST5)
st5_output$meanA[b] <- st5_output$meanG[b] <- mean(popPar_St5$gv_a)
st5_output$varG[b] <- varG(ST5)
st5_output$h2[b] <- varG(ST5) / varP(ST5) 

st5_het <- pullQtlGeno(ST5)
st5_het <- mean(rowMeans(1 - abs(st5_het - 1)))

st5_output$Fcoeff[b] <- 1 - (st5_het / st5_het_0)
st5_output$rEBV[b] <- cor(ST5@pheno, ST5@gv)

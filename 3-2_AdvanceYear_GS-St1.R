# Advance breeding program by 1 year
# Works backwards through pipeline to avoid overwriting data before selection.

# Year 7
# Release variety

# Year 6
ST5 <- setPheno(ST4,varE = varE_St1, reps = rep_S45, fixEff = fx_eff)
ST5@pheno <- (ST4@pheno + ST5@pheno) / 2

# Year 5
ST4 <- selectInd(ST3, n_S4)
ST4 <- setPheno(ST4, varE = varE_St1, reps = rep_S45, fixEff = fx_eff)

# Year 4
ST3 <- selectInd(ST2, n_S3)
ST3 <- setPheno(ST3, varE = varE_St1, reps = rep_S3, fixEff = fx_eff)

# Year 3
ST2 <- selectInd(ST1, n_S2, use = "ebv") # A: GEBV; A+D: GEGV                 
ST2 <- setPheno(ST2, varE = varE_St1, reps = rep_S2, fixEff = fx_eff)

# Year 2
ST1 <- selectInd(seedlings, n_S1)                         
ST1 <- setPheno(ST1, varE = varE_St1, reps = 1, fixEff = fx_eff)
                                                   



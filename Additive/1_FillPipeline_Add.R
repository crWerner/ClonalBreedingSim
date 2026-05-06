
# Set up initial yield trials with unique individuals

for (year in 1:6) {
  cat("Fill Pipeline year:", year, "of 6\n")
  
  if(year < 7){ # Year 1
    seedlings <- randCross(parents, n_crosses, nProgeny = n_famLines)
    seedlings <- setPheno(seedlings, varE = varE_seed, reps = 1)
  }
  
  if(year < 6){ # Year 2
    seedlings <- selectWithinFam(seedlings, 7)
    ST1 <- selectInd(seedlings, n_S1)
    ST1 <- setPheno(ST1, varE = varE_St1, rep = 1)
  }
  
  if(year < 5){ # Year 3
    ST2 <- selectInd(ST1, n_S2)
    ST2 <- setPheno(ST2, varE = varE_St1, reps = rep_S2)
  }
  
  if(year < 4){ #Y ear 4
    ST3 <- selectInd(ST2, n_S3)
    ST3 <- setPheno(ST3, varE = varE_St1, reps = rep_S3)
  }
  
  if(year<3){ # Year 5
    ST4 <- selectInd(ST3, n_S4)
    ST4 <- setPheno(ST4, varE = varE_St1, reps = rep_S45)
  }
  
  if(year < 2){ # Year 6
    ST5 <- setPheno(ST4, varE = varE_St1, reps = rep_S45)
    ST5@pheno = (ST4@pheno + ST5@pheno) / 2
  }
  
  if(year < 1){ # Year 7
  }
}

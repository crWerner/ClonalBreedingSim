
### EVERY YEAR, 20 OF 40 PARENTS ARE REPLACED BY
### 1. ALL ST3 - ST5 WHICH ARE NOT YET PART OF THE CROSSING POOL
### 2. ALL REMAINING SLOTS ARE FILLED WITH THE BEST PERFORMING S2 candidates 
###  --> regardless of whether or not they are better than the 30 parents discarded before.
###  --> reflects a general parental turnover of older parents.

#Drop Seedling parents
parents <- selectInd(pop = parents, nInd = (n_parents - n_new_parents))  # selects the best-performing current parents


#Add all S3-S5 lines which are not already included in nextParents from the years before
parent_cands <- c(ST3, ST4, ST5)
new_parents <- parent_cands[!parent_cands@id %in% parents@id]
rm(parent_cands)

if (nInd(new_parents) > n_new_parents) {
  new_parents <- selectInd(pop = new_parents, nInd = n_new_parents)
}

parents <- c(parents, new_parents)


# Fill all empty crossing block slots up to total number of parents with the best performing individuals from stage 2

if (nInd(parents) < n_parents) {
  st2_parents <- selectInd(pop = ST2, nInd = (n_parents - nInd(parents)))
  parents <- c(parents, st2_parents)
  rm(st2_parents)
}



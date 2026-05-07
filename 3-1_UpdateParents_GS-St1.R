
# For the sake of simplicity, new parental candidate only come from the 
# most recent Stage 1 population.

parent_cands <- c(parents, ST1)
parents <- selectInd(parent_cands, n_parents, use = "ebv")   # ebv = GEBV

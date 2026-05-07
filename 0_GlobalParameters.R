
# GENOME SIMULATION
species <- "GENERIC"   # see "?runMacs"
n_chr <- 10             # total number of chromosomes
n_qtl <- 200            # number of QTL per chromosome
n_snp <- 200            # number of SNP per chromosome (20k SNP array)

# TRAIT PARAMETERS
init_meanG <- 10        # genetic mean in the founder population       
init_varG <- 1          # genetic variance in the founder population

dominance_degree <- 0.9      # mean dominance degree. Set to 0 for additive trait.
dominance_degree_var <- 0.2  # variance of dominance degree

# TARGET PLOT-LEVEL HERITABILITIES 
H2_seed <- 0.1     # ... in the first seedlings
H2_St1 <- 0.3      # ... in the first clonal stage


### BREEDING PROGRAM PARAMETERS ###

# DURATION
brnIn <- 10     # number of "historic" (burn-in) breeding cycles
cycles <- 8    # number of "future" breeding cycles


# PARENTS AND CROSSINGS

# Parents 
n_founders <- 40     # number of founder individuals
n_parents <- 40      # total number of parents
n_new_parents <- 20  # new parents added each year to the crossing pool using phenotypic selection

# Crosses produced in the different scenarios
n_crosses <- 150          # number of crosses produced when no GS is used
n_crosses_GSseed <- 130   # number of crosses when GS is used 

# Lines produced per cross under the 4 different simulation scenarios
# When GS is used, the number of crosses and genotypes per cross is reduced
# to compensate for increased costs.
n_famLines <- 100          # number of progeny per cross : 150 crosses x 100 genotypes = 15.000
n_famLinesGS <- 90         # 150 x 90 = 13500
n_famLinesGSseed <- 80     # 130 x 80 = 10400
n_famLinesGSseed2 <- 40    # 130 x 40 x 2 = 10400; number of progeny per cross when 2 generation cycles per year are exploited

# Entries per stage (selection pressure)
n_S1 <- 1000
n_S2 <- 100
n_S3 <- 20
n_S4 <- 5

fam_maxSeed <- ceiling(n_S1 / n_crosses)  # max. number of advanced seedlings per family. Used "FillPipeline" for balanced family representation.
fam_maxSeed_GPCP <- 10

# Effective replications of yield trials to define heritabilities based on eVarSt1.
rep_S2 <- 2      # 2 reps in Stage 2
rep_S3 <- 4      # 4 reps in Stage 3
rep_S45 <- 6     # 6 reps in Stage 4 & 5


### FUNCTIONS

# Error variance required to obtain target Repeatability
var_e <- function(var_g, H2) {
  (var_g / H2) - var_g
}


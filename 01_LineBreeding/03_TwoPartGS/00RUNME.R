# Script name: Two-part wheat line breeding program
#
# Authors: Initially developed by Chris Gaynor; exanded/polished for this
# publication by Jon Bancic, Philip Greenspoon, Chris Gaynor, Gregor Gorjanc
#
# Date Created: 2023-01-23
#
# Uses two-part strategy with rapid cycling of parents in population improvement
# and conventional breeding for product development. Applies GS to advance
# individuals from DH to make PYT as well as in population improvement.

### empirically looks like these results might be plausible

#### standardize the burnin to be the same across all simulations, run for 40? reps,
##### change plots to ggplot and overlay into a single graph/plot (better viz)
#### have baseline of continuing pheno selection rather than GS

######## for later, once confident with simulations 
#### if we aren't retraining, we can do speed breeding/rapid cycling without phenotyping.
###### they only have 2 rapid cycles per year in the pop improvement

## musings from Dan
### maybe more important than model and retraining, would be to create an ideal population for long term gain
#### makes genomic selection and rapid cycling most effective
#### best pop for genomic prediction has intermediate num of segregating loci (relative to trait),
#### (accuracy = #indiv/#loci * heritability) large num loci means ridiculously large pop
#### F2/biparental population has 60 loci, but runs out of diversity fast
#### wild diversity population has 100s of millions of loci
#### design a population maybe MAGIC pop of 8 or so (which normally are only discussed for QTL mapping) --- 800ish effective loci, 1000 indivs, high heritability = 0.5 accuracy
#### time and cost are at play because num of gens/crossings to create MAGIC pop

#### once sims are good, get to year 0, build a MAGIC pop, then continue with GS

# ---- Clean environment and load packages ----
rm(list = ls())

#### probably load here the burn in pop

# install.packages(pkgs = "AlphaSimR")
#library(package = "AlphaSimR")
Xpackagers <- c('AlphaSimR','bWGR','parallel','foreach','doParallel',
                'reshape','ggplot2','gridExtra','lubridate','plyr',
                'ranger','Rcpp','keras','verification','rrBLUP',
                'reshape2','ScottKnott','viridis')
#install.packages(Xpackagers)
XXX <- lapply(Xpackagers, function(x){suppressMessages(require(x,quietly = TRUE, character.only = TRUE))})

# ---- Load global parameters ----
source(file = "GlobalParameters.R")
source(file="RunGSModels.R")
scenarioName = "LineGSTP_CS_noRetrain"
#scenarioName = "LineGSTP_BayesB_retrainNA"
#bayesB="emBB"

### retrain: TRUE (normal yearly), 5 (every 5 years), FALSE (none)
retrain = FALSE

# ---- Create list to store results from reps ----
results = list()
results_accPI = list()

for(REP in 1:nReps){
  cat("Working on REP:", REP,"\n")

  # ---- Create a data frame to track key parameters ----
  output = data.frame(year     = 1:nCycles,
                      rep      = rep(REP, nCycles),
                      scenario = "",
                      meanG = numeric(nCycles),
                      varG  = numeric(nCycles),
                      accSel  = numeric(nCycles))

  # ---- Create initial parents ----
  source(file = "CreateParents.R")

  # ---- Fill breeding pipeline with unique individuals from initial parents ----
  source(file = "FillPipeline.R")

  # ---- Simulate year effects ----
  P = runif(nCycles)

  # ---- Burn-in phase: Phenotypic selection program ----
  cat("--> Working on Burn-in \n")
  for(year in 1:nBurnin) {
    cat(" Working on burn-in year:",year,"\n")
    source(file = "UpdateParents.R") # Pick new parents
    source(file = "AdvanceYear.R")   # Advance yield trials by a year
    source(file = "StoreTrainPop.R") # Store training population
    # Report results
    output$meanG[year] = meanG(DH)
    output$varG[year]  = varG(DH)
  }

  # ---- Future phase: Genomic selection program ----
  cat("--> Working on Two-part line breeding program \n")
  # New parameters for population improvement
  nCyclesPI = 2    # Number of rapid cycles per year
  nParents  = 50   # Number of parents
  nCrossPI  = 100  # Number of crosses per cycle
  nF1PI = 100      # Number of F1-PI to advance to PD
  # Create a data frame to track selection accuracy in PI
  accPI = data.frame(accPI = numeric(nFuture*nCyclesPI))

  for(year in (nBurnin+1):(nBurnin+nFuture)) {
    cat(" Working on future year:",year,"\n")
    source(file = "RunGSModels.R")      # Run genomic model
    source(file = "AdvanceYear_GSTP.R")   # Advance yield trials by a year
    # Store training population
    if(retrain){source(file = "StoreTrainPop.R")} else if(is.numeric(retrain)){
      source(file = "StoreTrainPop_UpdateEvery5.R")
    } else{}  
    # Report results
    output$meanG[year] = meanG(DH)
    output$varG[year]  = varG(DH)
  }

  # Save results from current replicate
  results = append(results, list(output))
  results_accPI = append(results_accPI, list(accPI))
}

# Save results
saveRDS(results, file = paste0(scenarioName,".rds"))
saveRDS(results_accPI, file = paste0(scenarioName,"_accPI.rds"))

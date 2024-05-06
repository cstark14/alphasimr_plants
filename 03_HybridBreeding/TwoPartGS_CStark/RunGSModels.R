#### HAVE TO MAKE THIS INTO A FUNCTION. Otherwise, it won't be able to capture new DH

# Run genomic models
cat("  Running GS model\n")
i=1
if(exists("bayesB")){
  inputs <- vector(mode="list")
  modelsM <- vector(mode="list")
  modelsF <- vector(mode="list")
  ebvsM <- vector(mode="list")
  ebvsF <- vector(mode="list")
  
  #male_geno_MC <- data.frame(NULL);
  #male_pheno_MC <- data.frame(NULL);
  #female_geno_MC <- data.frame(NULL);
  #female_pheno_MC <- data.frame(NULL);
  
  maleTrainGen <- data.frame(i="Training",pullSnpGeno(MaleTrainPop,simParam = SP))
  maleYT1Gen <- data.frame(i="YT1",pullSnpGeno(MaleYT1,simParam = SP))
  maleDHGen <- data.frame(i="DH",pullSnpGeno(MaleDH,simParam = SP))
  malesPop <- c(MaleTrainPop,MaleYT1)
  malesGen <- pullSnpGeno(malesPop, simParam = SP)
  male_geno_MC <- rbind(data.frame(i=1,malesGen),male_geno_MC)
  male_pheno_MC <- rbind(data.frame(i=1,MaleTrainPop@pheno[,1]),male_pheno_MC)
  
  
  ### below would incorporate geno data for training and for new pop
  gsModelM = emBB(as.matrix(male_pheno_MC[,-1]),malesGen)
  
  if(anyNA(gsModelM$hat)){
    gsModelM$hat=rnorm(length(fit$hat))
  } 
  ### Need to figure out how to use full Genos but only get new pop data. 
  gsModelAllMales = gsModelM$hat
  ###Subset out the num in training pop?
  #pop@ebv <- as.matrix(fit$hat[pheno_MC$i == i],ncol=1)
  gsModelM <- as.matrix(gsModelAllMales[(nrow(male_pheno_MC)+1):length(gsModelAllMales)],ncol=1)
  
  femalesPop <- c(FemaleTrainPop,FemaleYT1)
  femalesGen <- pullSnpGeno(femalesPop, simParam = SP)
  female_geno_MC <- rbind(data.frame(i=1,femalesGen),female_geno_MC)
  female_pheno_MC <- rbind(data.frame(i=1,FemaleTrainPop@pheno[,1]),female_pheno_MC)

  gsModelF = emBB(as.matrix(female_pheno_MC[,-1]),as.matrix(female_geno_MC[,-1]))
  if(anyNA(gsModelF$hat)){
    gsModelF$hat=rnorm(length(fit$hat))
  }
  gsModelAllFemales = gsModelF$hat
  # pop@ebv <- as.matrix(fit$hat[pheno_MC$i == i],ncol=1)
  gsModelF <- as.matrix(gsModelAllFemales[(nrow(female_pheno_MC)+1):length(gsModelAllFemales)],ncol=1)
} else{
  
  gsModelM = RRBLUP(MaleTrainPop)
  gsModelF = RRBLUP(FemaleTrainPop)
}

### code that didn't work but keeping in a dummy function for future ideas
dummyGS <- function(){
  #maleGen = pullSnpGeno(MaleTrainPop, simParam = SP)
  #maleYT1Gen = pullSnpGeno(MaleYT1, simParam = SP)
  ### below would include geno data just for new pop... which doesn't seem right.
  #gsModelM = emBB(as.matrix(male_pheno_MC[,-1]),maleYT1Gen)
  
  #  male_pheno_MC <- rbind(data.frame(i=1,MaleTrainPop@pheno[,1]-mean(MaleTrainPop@pheno[,1])),male_pheno_MC)
  # female_pheno_MC <- rbind(data.frame(i=1,FemaleTrainPop@pheno[,1]-mean(FemaleTrainPop@pheno[,1])),female_pheno_MC)
  
  
  # got this code from the Features/Techniques section -> genomic models
  # Alternatively, one can use an external package
  # rrBLUP package is used for demonstration of GBLUP
  # y = pop@pheno[,1] 
  # G = A.mat(pullSnpGeno(pop, snpChip = 2)-1)
  # ans = mixed.solve(y = y, K = G)
  # cor(ans$u, pop@ebv) # results equivalent to RRBLUP 
  # pop@ebv = as.matrix(ans$u) # assign EBVs to population
}
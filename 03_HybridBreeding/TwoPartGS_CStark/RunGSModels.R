# Run genomic models
cat("  Running GS model\n")
i=1
if(exists("bayesB")){
  male_geno_MC <- data.frame(NULL);
  male_pheno_MC <- data.frame(NULL);
  
  female_geno_MC <- data.frame(NULL);
  female_pheno_MC <- data.frame(NULL);
  
  #maleGen = pullSnpGeno(MaleTrainPop, simParam = SP)
  #maleYT1Gen = pullSnpGeno(MaleYT1, simParam = SP)
  malesPop <- c(MaleTrainPop,MaleYT1)
  malesGen <- pullSnpGeno(malesPop, simParam = SP)
  male_geno_MC <- rbind(data.frame(i=1,maleGen),male_geno_MC)
  male_pheno_MC <- rbind(data.frame(i=1,MaleTrainPop@pheno[,1]-mean(MaleTrainPop@pheno[,1])),male_pheno_MC)
  
  ### below owuld include geno data just for new pop... which doesn't seem right.
  #gsModelM = emBB(as.matrix(male_pheno_MC[,-1]),maleYT1Gen)
  ### below would incorporate geno data for training and for new pop
  gsModelM = emBB(as.matrix(male_pheno_MC[,-1]),malesGen)
  
  # got this code from the Features/Techniques section -> genomic models
  #maleG = A.mat(pullSnpGeno(MaleTrainPop, snpChip = 2)-1)
  #ans = mixed.solve(y = y, K = G)
  #cor(ans$u, pop@ebv) # results equivalent to RRBLUP 
  #pop@ebv = as.matrix(ans$u) # assign EBVs to population
  
  if(anyNA(gsModelM$hat)){
    gsModelM$hat=rnorm(length(fit$hat))
  } 
  ### Need to figure out how to use full Genos but only get new pop data. 
  gsModelAllMales = gsModelM$hat
  ###Subset out the num in training pop?
  gsModelM <- as.matrix(gsModelAllMales[(nrow(male_pheno_MC)+1):length(gsModelAllMales)],ncol=1)
  
  femalesPop <- c(FemaleTrainPop,FemaleYT1)
  femalesGen <- pullSnpGeno(femalesPop, simParam = SP)
  female_geno_MC <- rbind(data.frame(i=1,femalesGen),female_geno_MC)
  female_pheno_MC <- rbind(data.frame(i=1,FemaleTrainPop@pheno[,1]-mean(FemaleTrainPop@pheno[,1])),female_pheno_MC)

  gsModelF = emBB(as.matrix(female_pheno_MC[,-1]),as.matrix(female_geno_MC[,-1]))
  if(anyNA(gsModelF$hat)){
    gsModelF$hat=rnorm(length(fit$hat))
  }
  ### Need to figure out how to use full Genos but only get new pop data. 
  gsModelAllFemales = gsModelF$hat
  ###Subset out the num in training pop?
  gsModelF <- as.matrix(gsModelAllFemales[(nrow(female_pheno_MC)+1):length(gsModelAllFemales)],ncol=1)
} else{
  
  gsModelM = RRBLUP(MaleTrainPop)
  gsModelF = RRBLUP(FemaleTrainPop)
}
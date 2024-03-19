# Run genomic models
cat("  Running GS model\n")
i=1
if(exists("bayesB")){
  male_geno_MC <- data.frame(NULL);
  male_pheno_MC <- data.frame(NULL);
  
  female_geno_MC <- data.frame(NULL);
  female_pheno_MC <- data.frame(NULL);
  
  maleGen = pullSnpGeno(MaleTrainPop, simParam = SP)
  male_geno_MC <- rbind(data.frame(i=1,maleGen),male_geno_MC)
  male_pheno_MC <- rbind(data.frame(i=1,MaleTrainPop@pheno[,1]-mean(MaleTrainPop@pheno[,1])),male_pheno_MC)
  
  gsModelM = emBB(as.matrix(male_pheno_MC[,-1]),as.matrix(male_geno_MC[,-1]))
  if(anyNA(gsModelM$hat)){
    gsModelM$hat=rnorm(length(fit$hat))
  } 
  gsModelM = as.matrix(gsModelM$hat[male_pheno_MC$i == i],ncol=1)
  
  femaleGen = pullSnpGeno(FemaleTrainPop, simParam = SP)
  female_geno_MC <- rbind(data.frame(i=1,femaleGen),female_geno_MC)
  female_pheno_MC <- rbind(data.frame(i=1,FemaleTrainPop@pheno[,1]-mean(FemaleTrainPop@pheno[,1])),female_pheno_MC)
  
  gsModelF = emBB(as.matrix(female_pheno_MC[,-1]),as.matrix(female_geno_MC[,-1]))
  if(anyNA(gsModelF$hat)){
    gsModelF$hat=rnorm(length(fit$hat))
  } 
  gsModelF = as.matrix(gsModelF$hat[female_pheno_MC$i == i],ncol=1)
} else{
  
  gsModelM = RRBLUP(MaleTrainPop)
  gsModelF = RRBLUP(FemaleTrainPop)
}
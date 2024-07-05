# Run genomic model

cat("  Running GS model\n")
#gsModel = RRBLUP(TrainPop)
runGS <- function(sex=NA,targetPop="YT1",useBayes=exists("bayesB")){
  if(is.na(sex)){
    if(useBayes==FALSE){
      gsModel = RRBLUP(TrainPop)
      return(gsModel)
    }
    TrainPheno <- data.frame(i="Training",TrainPop@pheno[,1])
    TrainGen <- data.frame(i="Training",pullSnpGeno(TrainPop,simParam = SP))
    
    if(targetPop=="YT1"){
      YT1Gen <- data.frame(i="YT1",pullSnpGeno(YT1,simParam = SP))
      Gen <- rbind(TrainGen,YT1Gen)
    } else if(targetPop=="DH"){
      DHGen <- data.frame(i="DH",pullSnpGeno(DH,simParam = SP))
      Gen <- rbind(TrainGen,DHGen)
    } else if(targetPop=="Parent"){
      ParentGen <- data.frame(i="Parent",pullSnpGeno(Parents,simParam = SP))
      Gen <- rbind(TrainGen,ParentGen)
    } else{
      cat(paste0("No target pop found for ",targetPop," . Please choose from YT1, DH, or Parent"))
      break
    }
    
    tries = 1
    # Repeat loop
    repeat {
      if (tries > 5) {
        break
      } 
      ### below would incorporate geno data for training and for new pop
      gsModel = emBB(as.matrix(TrainPheno[,-1]),as.matrix(Gen[,-1]))
      if(!is.na(gsModel$Ve)){break}
      # Increment x by 1
      tries = tries + 1
    }
    if(anyNA(gsModel$hat)){
      gsModel$hat=rnorm(length(gsModel$hat))
    } 
    
    ###Subset out the num in training pop?
    #pop@ebv <- as.matrix(fit$hat[pheno_MC$i == i],ncol=1)
    gsModel <- as.matrix(gsModel$hat[Gen$i==targetPop],ncol=1)
    return(gsModel)
  }else if(sex=="MALE"){
    if(useBayes==FALSE){
      gsModelM = RRBLUP(MaleTrainPop)
      return(gsModelM)
    }
    maleTrainPheno <- data.frame(i="Training",MaleTrainPop@pheno[,1])
    maleTrainGen <- data.frame(i="Training",pullSnpGeno(MaleTrainPop,simParam = SP))
    
    if(targetPop=="YT1"){
      maleYT1Gen <- data.frame(i="YT1",pullSnpGeno(MaleYT1,simParam = SP))
      malesGen <- rbind(maleTrainGen,maleYT1Gen)
    } else if(targetPop=="DH"){
      maleDHGen <- data.frame(i="DH",pullSnpGeno(MaleDH,simParam = SP))
      malesGen <- rbind(maleTrainGen,maleDHGen)
    } else if(targetPop=="Parent"){
      maleParentGen <- data.frame(i="Parent",pullSnpGeno(MaleParents,simParam = SP))
      malesGen <- rbind(maleTrainGen,maleParentGen)
    } else{
      cat(paste0("No target pop found for ",targetPop," . Please choose from YT1, DH, or Parent"))
      break
    }
    
    tries = 1
    # Repeat loop
    repeat {
      if (tries > 5) {
        break
      } 
      ### below would incorporate geno data for training and for new pop
      gsModelM = emBB(as.matrix(maleTrainPheno[,-1]),as.matrix(malesGen[,-1]))
      if(!is.na(gsModelM$Ve)){break}
      # Increment x by 1
      tries = tries + 1
    }
    if(anyNA(gsModelM$hat)){
      gsModelM$hat=rnorm(length(gsModelM$hat))
    } 
    
    ###Subset out the num in training pop?
    #pop@ebv <- as.matrix(fit$hat[pheno_MC$i == i],ncol=1)
    gsModelM <- as.matrix(gsModelM$hat[malesGen$i==targetPop],ncol=1)
    return(gsModelM)
  } else if(sex=="FEMALE"){
    if(useBayes==FALSE){
      gsModelM = RRBLUP(FemaleTrainPop)
      return(gsModelF)
    }
    femaleTrainPheno <- data.frame(i="Training",FemaleTrainPop@pheno[,1])
    femaleTrainGen <- data.frame(i="Training",pullSnpGeno(FemaleTrainPop,simParam = SP))
    
    if(targetPop=="YT1"){
      femaleYT1Gen <- data.frame(i="YT1",pullSnpGeno(FemaleYT1,simParam = SP))
      femalesGen <- rbind(femaleTrainGen,femaleYT1Gen)
    } else if(targetPop=="DH"){
      femaleDHGen <- data.frame(i="DH",pullSnpGeno(FemaleDH,simParam = SP))
      femalesGen <- rbind(femaleTrainGen,femaleDHGen)
    } else if(targetPop=="Parent"){
      femaleParentGen <- data.frame(i="Parent",pullSnpGeno(FemaleParents,simParam = SP))
      femalesGen <- rbind(femaleTrainGen,femaleParentGen)
    } else{
      print(paste0("No target pop found for ",targetPop," . Please choose from YT1, DH, or Parent"))
      break
    }
    
    tries = 1
    # Repeat loop
    repeat {
      if (tries > 5) {
        break
      } 
      ### below would incorporate geno data for training and for new pop
      gsModelF = emBB(as.matrix(femaleTrainPheno[,-1]),as.matrix(femalesGen[,-1]))
      if(!is.na(gsModelF$Ve)){break}
      # Increment x by 1
      tries = tries + 1
    }
    if(anyNA(gsModelF$hat)){
      gsModelF$hat=rnorm(length(gsModelF$hat))
    } 
    
    ###Subset out the num in training pop?
    #pop@ebv <- as.matrix(fit$hat[pheno_MC$i == i],ncol=1)
    gsModelF <- as.matrix(gsModelF$hat[femalesGen$i==targetPop],ncol=1)
    return(gsModelF)
  } else{
    cat(paste0("Sex (",sex,") did not match either MALE or FEMALE."))
  }
}

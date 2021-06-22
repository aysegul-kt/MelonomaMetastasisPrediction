prepareDataPool <- function(){
  Grp1 <-  ttest.results_detailed[ttest.results_detailed$dtype %in% "miRNA",]
  Grp1 <-  Grp1[Grp1$Diff001 %in% "1",]
  Grp1 <-  Grp1[Grp1$direction %in% c("L", "H"),]
  Grp1.L <-  Grp1[Grp1$direction %in% c( "L"),]
  Grp1.H <-  Grp1[Grp1$direction %in% c( "H"),]
  
  Grp2 <-  ttest.results_detailed[ttest.results_detailed$dtype %in% "mRNA",]
  Grp2 <-  Grp2[Grp2$Diff001 %in% "1",]
  Grp2 <-  Grp2[Grp2$direction %in% c("L", "H"),]
  Grp2.L <-  Grp2[Grp2$direction %in% c( "L"),]
  Grp2.H <-  Grp2[Grp2$direction %in% c( "H"),]
  
  Grp3 <-  GeneMethylationDetailedResults[GeneMethylationDetailedResults$dtype %in% "Methylation",]
  Grp3 <-  Grp3[Grp3$Diff001 %in% "1",]
  Grp3 <-  Grp3[Grp3$direction %in% c("L", "H"),]
  
  
  colnames(Grp3)<-colnames(Grp1)
 
 

  
  
  inputSignificanceSet <- rbind(Grp1,Grp2,Grp3) 
  
  #SelectedVar <- FilterSignificantVariables(inputSignificanceSet,TRUE)
  SelectedVar <- unique(inputSignificanceSet$probName)
  cat(paste("\n Variable Count",length(SelectedVar)))
  
  
  
  ### Cerate data from selected variables
  dataPool <- subset(finalDataSet,select = as.vector( SelectedVar ))
  cat("\n Step 1")
  dataPool <- dataPool[ , apply(dataPool, 2, function(x) !any(is.na(x)))]  # remove columns with NA
  
  Grp3 <- Grp3[Grp3$probName %in%  colnames(dataPool), ]
  Grp3.H <-  Grp3[Grp3$direction %in% c( "H"),]
  Grp3.L <-  Grp3[Grp3$direction %in% c( "L"),]
  cat("\n Step 2- grp 3 modified")
  
  dataPool$class <- finalDataSet$class
  
  dataPool.class1<-dataPool[dataPool$class %in% "Metastatic", ]  
  dataPool.class2<-dataPool[dataPool$class %in% "Primary Tumor", ]
  cat("\n Step 3")
  dataPool.class1$class<-rep(c(1),nrow(dataPool.class1))
  dataPool.class2$class<-rep(c(0),nrow(dataPool.class2))
  cat("\n Step 4")
  
  # 1: metastastic 0: Primary tumor
  dataPool.class1$class<-factor(dataPool.class1$class, levels = c(0, 1)) 
  dataPool.class2$class<-factor(dataPool.class2$class, levels = c(0, 1)) 
  dataPool <- rbind(dataPool.class1,dataPool.class2)
  dataPool <- dataPool[sample(1:nrow(dataPool), nrow(dataPool), replace=FALSE),] 
  
  cat(paste("\n Data Pool Size:",ncol(dataPool)))
  
  response <- NULL
  response$data<- dataPool
  
  response$mRNAFeatures <- Grp2
  response$mRNAFeatures.H <- Grp2.H
  response$mRNAFeatures.L <- Grp2.L
  
  response$miRNAFeatures <- Grp1
  response$miRNAFeatures.H <- Grp1.H
  response$miRNAFeatures.L <- Grp1.L
  
  response$methylationFeatures  <- Grp3
  response$methylationFeatures.H  <- Grp3.H
  response$methylationFeatures.L  <- Grp3.L
  
  response
  
}
getSubSampleOfDataPool <- function(dataPool){
  response <- NULL
  datasubsetclass1 <- dataPool[dataPool$class %in% 1, ]  
  datasubsetclass2  <- dataPool[dataPool$class %in% 0, ]  
  
  s <- sample(1:nrow(datasubsetclass1), nrow(datasubsetclass2),replace=FALSE)
  class1 <- datasubsetclass1[s,] 
  class1Remain <- datasubsetclass1[-s, ] 
  
  # class1 <-  class1[sample(1:nrow(class1), nrow(class1), replace=FALSE),]
  trainInd <- round(nrow(class1)*0.8) 
  # create test and train sets for class 1
  class1.train_set <- class1[1:trainInd,]
  class1.test_set <- class1[-(1:trainInd),]
  
  
  trainInd <- round(nrow(class1Remain)*0.8) 
  
  class1Remain.train_set <- class1Remain[1:trainInd,]
  class1Remain.test_set <- class1Remain[-(1:trainInd),]
  
  
  #reorder the class 2 and split from %80
  class2 <- datasubsetclass2[sample(1:nrow(datasubsetclass2), nrow(datasubsetclass2), replace=FALSE),]
  
  trainInd <- round(nrow(class2)*0.8) 
  # create test and train sets for class 2
  class2.train_set<-class2[1:trainInd,]
  class2.test_set <- class2[-(1:trainInd),]
  
  #bind traning sets for class 1 and class 2
  data.training <- rbind(class1.train_set,class2.train_set)
  #re order the sets
  data.training <- data.training[sample(1:nrow(data.training), nrow(data.training), replace=FALSE),]
  
  #Commet out for Smote
   set <- rbind(class1Remain.train_set,data.training)
   classIndex <- ncol(set)
   
   generatedData <- SMOTE(set[,-c(classIndex)],set[,classIndex], K=5)
   generatedData <-as.data.frame( generatedData$data )
  
  #bind traning and test set for class 1 and class 2
  data.test <- rbind(class2.test_set,class1.test_set)
  # re order the test set
  data.test <- data.test[sample(1:nrow(data.test), nrow(data.test), replace=FALSE),]
  
  response$data.training <- data.training
  response$data.test <- data.test
  response$remain.train_set <- class1Remain.train_set
  response$remain.test_set <- class1Remain.test_set
  response$smote.training <- generatedData
 
  response
  
}
getPCAPredictions <- function (data.training,data.test, isAll) {
  ################ PCA ############################
  # principal component Analaysis 
  
  data.training.pca<- prcomp(data.training[,-c(ncol(data.training))])
  
  # Comment out for PCA Logs and figures..
  #compute standard deviation of each principal component
  std_dev <- data.training.pca$sdev
  #compute variance
  pr_var <- std_dev^2
  #check variance of first 100 components
  # pr_var[1:100]
  #proportion of variance explained
  prop_varex <- pr_var/sum(pr_var)
  sumcomp <-0
  for (comp in 1: length( prop_varex)) {
    sumcomp <- sumcomp + prop_varex[comp]
    if (sumcomp < 0.98) { XX <- comp }
  }
  
  # Prepare for prediction
  data.test.pca <- predict(data.training.pca, newdata = data.test)
  
  
  ############# MODEL RUN ##############
  
  cat( "\n PCA is  Ready")
  model.data.train <- data.frame(class = data.training$class , data.training.pca$x)
  #transform test into PCA
  
  # re order dataset
  model.data.train  <- model.data.train[sample(1:nrow(model.data.train), nrow(model.data.train), replace=FALSE),]
  data.test  <- data.test[sample(1:nrow(data.test), nrow(data.test), replace=FALSE),]
  
  # under sampling for metastastic to generate new set 
  class1<-model.data.train[model.data.train$class %in% 1, ]  
  class2<-model.data.train[model.data.train$class %in% 0, ]
 
  model.data.train<-rbind(class1,class2)
  
  class1<-data.test[data.test$class %in% 1, ]  
  class2<-data.test[data.test$class %in% 0, ]
 
  data.test<-rbind(class1,class2)
  test.data_outcomes <- data.test$class
  
  test.data <- predict(data.training.pca, newdata = data.test)
  test.data <- as.data.frame(test.data)
  test.data_outcomes<-data.test$class
  if (!isAll ) XX <- min(300, XX)
  
  cat(paste("\n Component Count:", XX,"\n"))
  response <-NULL
  response$trainning <- model.data.train[,1:XX]
  response$test <- test.data[,1:XX]
  response$test_outcomes <-test.data_outcomes
  response$importance <-NULL
  response
  
  
  
}


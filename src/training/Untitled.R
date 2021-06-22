getSubSampleOfDataPool <- function(datasubsetclass1,datasubsetclass2){
  response <- NULL

  class1 <- datasubsetclass1[sample(1:nrow(datasubsetclass1), nrow(datasubsetclass2), replace=FALSE),] 
  # class1 <-  class1[sample(1:nrow(class1), nrow(class1), replace=FALSE),]
  trainInd <- round(nrow(class1)*0.8) 
  # create test and train sets for class 1
  class1.train_set <- class1[1:trainInd,]
  class1.test_set <- class1[-(1:trainInd),]
  
  #reorder the class 2 and split from %80
  class2 <- datasubsetclass2[sample(1:nrow(datasubsetclass2), nrow(datasubsetclass2), replace=FALSE),]
  
  trainInd <- round(nrow(class2)*0.8) 
  # create test and train sets for class 2
  class2.train_set<-class2[1:trainInd,]
  class2.test_set <- class2[-(1:trainInd),]
  
  #bind traning sets for class 1 and class 2
  data.training <- rbind(class1.train_set,class2.train_set)
  #re order the sets
  data.training <- data.training_f[sample(1:nrow(data.training), nrow(data.training), replace=FALSE),]
  
  #Commet out for Smote
  # generatedData <- SMOTE(data.training[,-c(58978)],data.training[,58978], K=3)
  # data.training <-as.data.frame( generatedData$data )
  
  #bind traning and test set for class 1 and class 2
  data.test <- rbind(class2.test_set,class1.test_set)
  # re order the test set
  data.test <- data.test[sample(1:nrow(data.test), nrow(data.test), replace=FALSE),]
  
  response$data.training <- data.training
  response$data.test <- data.test
  
  
  response
  
}


getPCAPredictions <- function (data.training,data.test) {
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
  
  cat(paste("Try:",i, "Set:", t, "Metod:", k, j,  "\n"))
  model.data.train <- data.frame(Class = data.training$class , data.training.pca$x)
  #transform test into PCA
  
  # re order dataset
  model.data.train  <- model.data.train[sample(1:nrow(model.data.train), nrow(model.data.train), replace=FALSE),]
  data.test  <- data.test[sample(1:nrow(data.test), nrow(data.test), replace=FALSE),]
  
  # under sampling for metastastic to generate new set 
  class1<-model.data.train[model.data.train$Class %in% 1, ]  
  class2<-model.data.train[model.data.train$Class %in% 0, ]
  class1_selected <- class1[sample(1:nrow(class1), nrow(class1), replace=FALSE),]
  model.data.train<-rbind(class1_selected,class2)
  
  class1<-data.test[data.test$class %in% 1, ]  
  class2<-data.test[data.test$class %in% 0, ]
  class1_selected <- class1[sample(1:nrow(class1), nrow(class1), replace=FALSE),]
  data.test<-rbind(class1_selected,class2)
  test.data_outcomes <- data.test$class
  
  test.data <- predict(data.training.pca, newdata = data.test)
  test.data <- as.data.frame(test.data)
  test.data_outcomes<-data.test$class
  
  XX <- min(pc_count[j], XX)
  
  cat(paste("\n Component Count:", XX,"\n"))
  response <-NULL
  response$model.data.train <- model.data.train[,1:XX]
  response$test.data <- test.data[,1:XX]
  response$test.data_outcomes <-test.data_outcomes
  
  response
  
  
  
}

########################Funtions###############
fxDataPool<- function(ncount) {
  dataSample <-  getSubSampleOfDataPool(pool$data)
  dataSample
}

fxRfFilter <- function(training.data,test.data,varlist,usePrev,data){
  require(ranger)
  
  if(!usePrev){
    rg.input<- subset(training.data, select = varlist)
    rg.input$class <- training.data$class
    
    colnames.new <- lapply(colnames(rg.input),  function(x) {str_replace_all(x,"-","_x_")})
    colnames(rg.input) <- unlist(colnames.new)
    
    rf.data <- ranger(class ~ ., data = rg.input, importance = 'permutation')
    rg.importance <- importance_pvalues(rf.data, method = "altmann", formula = class ~ ., data = rg.input)
    
    imp.dt <- as.data.frame(rg.importance) 
    
    imp.dt$columnName <-rownames(imp.dt)
  }else {
    imp.dt <-data
  }
 

  
  fltr <- imp.dt[c(imp.dt$importance>mean(imp.dt$importance)) , ]
  cat("\nImportance Cut Of:")
  cat(mean(imp.dt$importance))
  fltr <- fltr[c(fltr$pvalue<"0.01") , ]
  
  selection <- unlist(lapply(fltr$columnName,  function(x) {str_replace_all(x,"_x_","-")}))
  
  selection <- selection[selection %in% colnames(training.data)]
  
  cat(paste("\n Component Count:", length(selection),"\n"))
  response <-NULL
  response$trainning <- subset(training.data, select = selection)
  
  response$trainning$class <- as.factor(training.data$class)### To do bu kısım miRNA run da sorun çıkarır 
  colnames(response$trainning) <- unlist(lapply(colnames(response$trainning),  function(x) {str_replace_all(x,"-","_x_")}))
  
  response$test <- subset(test.data, select = selection)
  colnames(response$test) <- unlist(lapply(colnames(response$test),  function(x) {str_replace_all(x,"-","_x_")}))
  
  response$test_outcomes <-test.data$class
  response$importance <- imp.dt
  response
  
}

fxGetModelDataInputs<- function (dataSample,mode,variableList,runtypes,usePrev,data,preruuid) {
  
  temp.traningSet <-  subset(dataSample$data.training, select = as.vector( variableList) )
  temp.traningSet$class <- dataSample$data.training$class
  
  cat(paste("columns: ", ncol(temp.traningSet)))
  temp.testSet <- subset(dataSample$data.test, select = as.vector( variableList) ) 
  temp.testSet$class <- dataSample$data.test$class
  
  tempRemoved.traningSet <- subset(dataSample$remain.train_set, select = as.vector(variableList))
  tempRemoved.traningSet$class <- dataSample$remain.train_set$class
  
  tempRemoved.testSet <- subset(dataSample$remain.test_set, select = as.vector(variableList ))
  tempRemoved.testSet$class <- dataSample$remain.test_set$class
  
  tempSmote.trainingSet <- subset(dataSample$smote.training, select = as.vector(variableList ))
  tempSmote.trainingSet$class <- dataSample$smote.training$class
  
  ModelInputs <-NULL
  
  if("us" %in% runtypes){
    # Under Sample
    if(mode =="PCA") modelData <- getPCAPredictions(temp.traningSet,temp.testSet, FALSE)
    else if(mode=="rf")  
    {
      importanceData <-NULL
      
      if(usePrev) {importanceData <- data[data$run_uuid %in% paste(preruuid, "-us"),]}
      
      modelData <- fxRfFilter(temp.traningSet,temp.testSet,variableList,usePrev = usePrev,data = importanceData)
      
      # rfResults$test$class <- rfResults$test_outcomes
      # modelData <- getPCAPredictions(rfResults$trainning,rfResults$test, TRUE)
      # modelData$importance <- rfResults$importance
    
    } else if( mode=="vars") {
      modelData <-NULL
     
      modelData$trainning <- subset(temp.traningSet, select = variableList)
      modelData$trainning$class <- as.factor(temp.traningSet$class)### To do bu kısım miRNA run da sorun çıkarır 
      
      colnames(modelData$trainning) <- unlist(lapply(colnames(modelData$trainning),  function(x) {str_replace_all(x,"-","_x_")}))
      
      modelData$test <- subset(temp.testSet, select = variableList)
      colnames(modelData$test) <- unlist(lapply(colnames(modelData$test),  function(x) {str_replace_all(x,"-","_x_")}))
      
      modelData$test_outcomes <-temp.testSet$class
      
    }
    ModelInputs$us <- modelData
    
  }
    
  if("all" %in% runtypes){
  
     # All Sample
    if(mode =="PCA")  modelData <- getPCAPredictions(rbind(temp.traningSet,tempRemoved.traningSet ) , rbind(temp.testSet,tempRemoved.testSet), FALSE)
    else if(mode=="rf")
    {
      importanceData <-NULL
      
      if(usePrev) {importanceData <- data[data$run_uuid %in% paste(preruuid, "-all"),]}
      modelData <- fxRfFilter(rbind(temp.traningSet,tempRemoved.traningSet ) , rbind(temp.testSet,tempRemoved.testSet),variableList,usePrev = usePrev,data = importanceData)
      #rfResults$test$class <- rfResults$test_outcomes
      #modelData <- getPCAPredictions(rfResults$trainning,rfResults$test, TRUE)
      #modelData$importance <- rfResults$importance
      #modelData <-rfResults
    }
    ModelInputs$all <- modelData
  }
  
  
  if("sm" %in% runtypes){
      
    # Smote 
    if(mode =="PCA")  modelData <- getPCAPredictions(tempSmote.trainingSet, temp.testSet, FALSE)
    else if(mode=="rf")   {
      importanceData <-NULL
      
      if(usePrev) {importanceData <- data[data$run_uuid %in% paste(preruuid, "-sm"),]}
      
      modelData <- fxRfFilter(tempSmote.trainingSet, temp.testSet,variableList,usePrev = usePrev,data = importanceData)
      # modelData <- rfResults
      #rfResults$test$class <- rfResults$test_outcomes
     # modelData <- getPCAPredictions(rfResults$trainning,rfResults$test,TRUE)
     # modelData$importance <- rfResults$importance
      #modelData <-rfResults
    } else if( mode=="vars") {
      modelData <-NULL
      
      modelData$trainning <- subset(tempSmote.trainingSet, select = variableList)
      modelData$trainning$class <- as.factor(tempSmote.trainingSet$class)### To do bu kısım miRNA run da sorun çıkarır 
      
      colnames(modelData$trainning) <- unlist(lapply(colnames(modelData$trainning),  function(x) {str_replace_all(x,"-","_x_")}))
      
      modelData$test <- subset(temp.testSet, select = variableList)
      colnames(modelData$test) <- unlist(lapply(colnames(modelData$test),  function(x) {str_replace_all(x,"-","_x_")}))
      
      modelData$test_outcomes <-temp.testSet$class
      
    }
    ModelInputs$sm <- modelData
  }
  return (ModelInputs)
  
}

fxPredict<- function (modeltype,run, modelInputs, Definition, variableList ,run_uuid,mode,runtypes) {
  
  Result_Matrix <- NULL
  importanceVals <-NULL
  
  model_list<-c(modeltype) # rbf eklenecek... "
  results <-NULL
  bestTunes <- NULL
  modelInfo <- NULL
  importanceVals <- NULL
  #UnderSample
  #PCA
  if ("us"  %in% runtypes)
  {
    cat(paste("\n run no:",run,"\n UnderSample run Starting"))
    
    modelData <- modelInputs$us
    
    cat(paste( "\n sample size",nrow(modelData$trainning) ))
    
    
    # Prediction
    exact <- system.time(Sonuc_Matrix_try <- trainModels( modelData$trainning, modelData$test, model_list,modelData$test_outcomes ))
    cat("\n")
    print(exact)
    cat("\n")
    # Add Logs to console
    print( ( Sonuc_Matrix_try$results[,c(1,2,12,17)]))
    
    
    # collect results
    Sonuc_Matrix_try$results <-
      fxAppendDefinitionsDetailes(
        Sonuc_Matrix_try$results,
        definition = paste(Definition, " : H+L: Undersample"),
        run_uuid = paste(run_uuid,"-us"),
        pcCount = ncol(modelData$test),
        run=run
      )
    results <-  Sonuc_Matrix_try$results
    
    Sonuc_Matrix_try$model$run <- rep(run,nrow(Sonuc_Matrix_try$model))
    Sonuc_Matrix_try$model$run_uuid <- rep(paste(run_uuid,"-us"),nrow(Sonuc_Matrix_try$model))
    modelInfo<- Sonuc_Matrix_try$model
    
    Sonuc_Matrix_try$bestTunes$run <- rep(run,nrow(Sonuc_Matrix_try$bestTunes))
    Sonuc_Matrix_try$bestTunes$run_uuid <- rep(paste(run_uuid,"-us"),nrow(Sonuc_Matrix_try$bestTunes))
    bestTunes<- Sonuc_Matrix_try$bestTunes
    
    if(mode=="rf")
    {
      modelData$importance$run_uuid <- rep(paste(run_uuid,"-us"),nrow(modelData$importance))
      modelData$importance$run <- rep(run,nrow(modelData$importance))
      importanceVals <-  modelData$importance
    }
  }
  
  #####################ALL 
  
  #PCA
  if("all" %in% runtypes)
  {
    
    # 
    cat(paste("\n run no:",run,"\n All Sample run Starting"))
    
    modelData <- modelInputs$all
    cat(paste( " \n sample size",nrow(modelData$trainning) ))
    
    #Predictions
    exact <- system.time(Sonuc_Matrix_try <- trainModels( modelData$trainning, modelData$test, model_list,modelData$test_outcomes ))
    cat("\n")
    print(exact)
    
    # # Add Logs to console
    cat("\n")
    print( ( Sonuc_Matrix_try$results[,c(1,2,12,17)]))
    
    #  # collect results
    Sonuc_Matrix_try$results <-
      fxAppendDefinitionsDetailes(
        Sonuc_Matrix_try$results,
        definition = paste(Definition, " : H+L: All"),
        run_uuid = paste(run_uuid,"-all"),
        pcCount = ncol(modelData$test),
        run=run
      )
    
    results <-  rbind (results,Sonuc_Matrix_try$results)
    
    Sonuc_Matrix_try$model$run <- rep(run,nrow(Sonuc_Matrix_try$model))
    Sonuc_Matrix_try$model$run_uuid <- rep(paste(run_uuid,"-all"),nrow(Sonuc_Matrix_try$model))
    modelInfo<- rbind (modelInfo,Sonuc_Matrix_try$model)
    
    Sonuc_Matrix_try$bestTunes$run <- rep(run,nrow(Sonuc_Matrix_try$bestTunes))
    Sonuc_Matrix_try$bestTunes$run_uuid <- rep(paste(run_uuid,"-all"),nrow(Sonuc_Matrix_try$bestTunes))
    bestTunes<- rbind (bestTunes,Sonuc_Matrix_try$bestTunes)
    
    if(mode=="rf")
    {
      modelData$importance$run_uuid <- rep(paste(run_uuid,"-all"),nrow(modelData$importance))
      modelData$importance$run <- rep(run,nrow(modelData$importance))
      importanceVals <-  rbind(importanceVals,modelData$importance)
    } 
  }
  
  
  if("sm" %in% runtypes){
    ########## SMOTE
    # SMOTE training
    # PCA
    cat(paste("\n run no:",run,"\n SMOTE run Starting"))
    
    modelData <- modelInputs$sm
    
    
    #Predictions
    exact <- system.time(Sonuc_Matrix_try <- trainModels( modelData$trainning, modelData$test, model_list,modelData$test_outcomes ))
    cat("\n")
    print(exact)
    
    # Log For results
    cat(paste( " \n sample size",nrow(modelData$trainning) ))
    print( ( Sonuc_Matrix_try$results[,c(1,2,12,17)]))
    
    # collect results
    Sonuc_Matrix_try$results <-
      fxAppendDefinitionsDetailes(
        Sonuc_Matrix_try$results,
        definition = paste(Definition, " : H+L: SMOTE"),
        run_uuid = paste(run_uuid,"-sm"),
        pcCount = ncol(modelData$test),
        run=run
      )
    results <-  rbind (results,Sonuc_Matrix_try$results)
    
    Sonuc_Matrix_try$model$run <- rep(run,nrow(Sonuc_Matrix_try$model))
    Sonuc_Matrix_try$model$run_uuid <- rep(paste(run_uuid,"-sm"),nrow(Sonuc_Matrix_try$model))
    modelInfo<- rbind (modelInfo,Sonuc_Matrix_try$model)
    
    Sonuc_Matrix_try$bestTunes$run <- rep(run,nrow(Sonuc_Matrix_try$bestTunes))
    Sonuc_Matrix_try$bestTunes$run_uuid <- rep(paste(run_uuid,"-sm"),nrow(Sonuc_Matrix_try$bestTunes))
    bestTunes<- rbind (bestTunes,Sonuc_Matrix_try$bestTunes)
    
    
    if(mode=="rf")
    {
      modelData$importance$run_uuid <- rep(paste(run_uuid,"-sm"),nrow(modelData$importance))
      modelData$importance$run <- rep(run,nrow(modelData$importance))
      importanceVals <-  rbind(importanceVals,modelData$importance)
    }
    
  }
  
  
  
  
  
  
  
  
  
  modelInfo$modelType <-  rep(modeltype,nrow(modelInfo)) 
  
  if(mode=="rf") {importanceVals$modelType <-  rep("common",nrow(importanceVals)) }
  bestTunes$modelType <-  rep(modeltype,nrow(bestTunes)) 
  
  Result_Matrix$modelInfo <- modelInfo
  Result_Matrix$results <- results 
  if(mode=="rf") { Result_Matrix$importance <-importanceVals}
  Result_Matrix$bestTunes <- bestTunes 
  Result_Matrix
  
}

fxAppendDefinitionsDetailes <- function(result,definition,run_uuid,pcCount,run){
  l<-length(length(result[1]))
  result$definition.type  <- rep(c(definition), l)
  result$definition.p <-  rep("0.001", l)
  result$definition.pcCount <- rep(pcCount, l)
  result$trialCode <- rep(run_uuid,l ) 
  result$dataset <- run
  result
}

calculateCutOf <- function(importance) {
  
  fltr <- importance[c(importance$importance>mean(importance$importance)) , ]
  cat("\nImportance Cut Of:")
  cat(mean(importance$importance))
  fltr <- fltr[c(fltr$pvalue<"0.01") , ]
  
  selection <- unlist(lapply(fltr$columnName,  function(x) {str_replace_all(x,"_x_","-")}))
  
  selection
}

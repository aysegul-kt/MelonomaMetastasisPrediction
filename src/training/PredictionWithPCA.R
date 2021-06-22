source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/training/PrepareDataPool.R')
library(smotefamily)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(ranger)

numCores <- detectCores()
numCores

pool  <- prepareDataPool()

View(pool)

system.time(dataSampleList <- lapply(rep(1,10), fxDataPool))

save(dataSampleList,file="data/processed/dataSampleList_Backup_1504.RData")

Result_Matrix_April <- NULL
Result_Matrix_April_Rf <-NULL
##############Run uuid#############
run_uuid <- UUIDgenerate(use.time = NA, n=1L)
#"9b7965fc-1cb1-4b31-8b90-c03b72e2e5ab"
#  how to access: dataSampleList[[1]][["data.training"]]
# run for MiRNA
# svmLinear,ranger, svmPoly,nnet,naive_bayes, svmRadialGrid,adaBoost
############Model List##############
models <- c( "svmLinear","ranger", "svmPoly","nnet","svmRadialGrid","svmRadialGrid3","adaBoost")
#models <- c("svmLinear")

#,as.vector(pool$mRNAFeatures$probName),as.vector(pool$methylationFeatures$probName)
####################Predictions################
## 
response <- NULL
run_uuid 

#variables <- c(as.vector(pool$miRNAFeatures$probName), as.vector(pool$methylationFeatures$probName))
lvariables <- rbind(pool$miRNAFeatures)
setorder(lvariables,P_Val) 

lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
lvariablesList <- as.vector(lvariables$probName)  
View(lvariablesList)


#run_uuid <- "a537332c-5a59-4b80-912d-190794996708"

#lvariablesList <-  selection

for (runn in c(1:10)) {

 #data <- Result_Matrix_April_Rf$importance[Result_Matrix_April_Rf$importance$run %in% runn,]  

modelInpts <-  fxGetModelDataInputs (dataSample=dataSampleList[[runn]],mode="rf",variableList=lvariablesList,runtypes=c("sm","us"),
                                     usePrev=FALSE,data=data,preruuid=NULL)
  
response <-
    lapply(
      models,
      fxPredict,
      run=runn,
      modelInputs=modelInpts ,
      Definition = "PCA_mRNA-Methy.L_0.001_base_500_Var", 
      variableList = lvariablesList,
      run_uuid = run_uuid,
      mode="rf", # PCA or rf, vars
      runtypes=c("sm","us")
    )

 
     Result_Matrix_April_Rf$importance <- rbind(Result_Matrix_April_Rf$importance,response[[1]]$importance)
  
  for (i  in 1:length(response)) {

    Result_Matrix_April_Rf$modelTrainResultsInfo <- rbindlist(list(Result_Matrix_April_Rf$modelTrainResultsInfo,response[[i]]$modelInfo), use.names=TRUE, fill=TRUE)
    Result_Matrix_April_Rf$result <-  rbind(Result_Matrix_April_Rf$result,response[[i]]$results)
    Result_Matrix_April_Rf$bestTune <- rbindlist(list(Result_Matrix_April_Rf$bestTune,response[[i]]$bestTunes), use.names=TRUE, fill=TRUE)
  }
    
    
}

######################## REPORT RESULTS #####################

View(Result_Matrix_April_Rf$result)

save(Result_Matrix_April_Rf,file="data/processed/Result_Matrix_May_Rf_2405.RData")
save(ResultsPCA,file="data/processed/ResultsPCA_Rf_2405.RData")


output <- CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"0c4aba8a-6506-4ce0-9a94-261934682414 -us","miRNA Only - Undersample(Hibrit)"  )
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"0c4aba8a-6506-4ce0-9a94-261934682414 -all"," miRNA Only- All (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"0c4aba8a-6506-4ce0-9a94-261934682414 -sm","miRNA Only- Smote (Hibrit)"  )) 

output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"9c11e543-d122-4167-becf-c35673eba0fc -us","miRNA & mRNA - Undersampl (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"9c11e543-d122-4167-becf-c35673eba0fc -all","miRNA & mRNA, - All (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"9c11e543-d122-4167-becf-c35673eba0fc -sm","miRNA & mRNA, - Smote (Hibrit)"  )) 

output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"a537332c-5a59-4b80-912d-190794996708 -us","miRNA , mRNA & Methy. - Undersample (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"a537332c-5a59-4b80-912d-190794996708 -all","miRNA , mRNA & Methy.- All (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"a537332c-5a59-4b80-912d-190794996708 -sm","miRNA , mRNA & Methy  Smote (Hibrit)"  )) 

output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"9f9e57a0-b82b-4ad2-a585-087f9e607dc9 -us","miRNA , mRNA & Methy L. Only- Undersample (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"9f9e57a0-b82b-4ad2-a585-087f9e607dc9 -sm","miRNA , mRNA & Methy - L  Smote (Hibrit)"  )) 

output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"8a1bfed7-31e7-4b0c-8f81-3be8b2fe8cb3 -us","miRNA , mRNA & Methy H. Only- Undersample (Hibrit)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"8a1bfed7-31e7-4b0c-8f81-3be8b2fe8cb3 -sm","miRNA , mRNA & Methy - H  Smote (Hibrit)"  )) 



#ResultsPCA <- NULL
ResultsPCA$result <- PCA_logs 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"b38ec272-9131-4838-b340-8e5e421d8562 -us","miRNA Only- Undersample (PCA)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"b38ec272-9131-4838-b340-8e5e421d8562 -all"," miRNA Only- All (PCA)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"b38ec272-9131-4838-b340-8e5e421d8562 -sm","miRNA, Only- Smote (PCA)"  )) 

output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"ab4cc487-2d45-46eb-9784-4449a731be27 -us","miRNA & mRNA - Undersampl (PCA)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"ab4cc487-2d45-46eb-9784-4449a731be27 -all","miRNA & mRNA, - All (PCA)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"ab4cc487-2d45-46eb-9784-4449a731be27 -sm","miRNA & mRNA, - Smote (PCA)"  )) 

output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"69fc2857-0534-4445-a7db-1dfba6a3a6d4 -us","miRNA , mRNA & Methy. Only- Undersample (PCA)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"69fc2857-0534-4445-a7db-1dfba6a3a6d4 -all","miRNA , mRNA & Methy.- All (PCA)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(ResultsPCA,"69fc2857-0534-4445-a7db-1dfba6a3a6d4 -sm","miRNA , mRNA & Methy  Smote (PCA)"  )) 

## after Tunning 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"372fdfbc-28e7-49c5-a9b3-0f691bb5c7b2 -us","miRNA , mRNA   Undersample (PCA revised) us"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"fd5685f1-f3e9-4fb1-bbc3-3cc6a0a687f0 -us","miRNA , mRNA   Undersample (PCA revised) us"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"fd5685f1-f3e9-4fb1-bbc3-3cc6a0a687f0 -sm","miRNA , mRNA   Undersample (PCA revised) sm"  )) 



## RF Tunning poly.... bad

output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"e8f88c0a-5d93-45a6-9779-e411e840a4ac -us","miRNA , mRNA   Undersample - US (RF revised)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"e8f88c0a-5d93-45a6-9779-e411e840a4ac -sm","miRNA , mRNA   Undersample- SM  (RF revises)"  )) 

# PCA smote retiry with new grid--- Bad
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"bdfe95f2-ed48-4b2a-9143-d422f9c06cbb -us","miRNA , mRNA   Undersample (PCA revised new grid) us"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"bdfe95f2-ed48-4b2a-9143-d422f9c06cbb -sm","miRNA , mRNA   Undersample (PCA revises new grid) sm"  )) 

# 2000 var with different tunnining- Nat Included>>> Modified all **** selected item u sette
# '1318c3f5-91e5-469b-b230-bb3e8de6e836'
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"1318c3f5-91e5-469b-b230-bb3e8de6e836 -us","miRNA , mRNA & Methy  PCA  Undersample (2000 var revised new grid) us"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"1318c3f5-91e5-469b-b230-bb3e8de6e836 -sm","miRNA , mRNA  & Methy PCA Smote (2000 var revises new grid) sm"  )) 

# RF tune 2000 vars - Nat Included  ???? 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"54507440-70e4-4b1b-a5e5-0c5cd890dd9a -us","miRNA , mRNA   Undersample (2000 var revised new grid-2) us"  )) 

# PCA - L methy 2000 vars -
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"3ae16ed8-43b1-4d0d-b0ca-08a0ec0f7341 -us","miRNA , mRNA, Methy L PCA Undersample (2000 var revised new grid)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"3ae16ed8-43b1-4d0d-b0ca-08a0ec0f7341 -sm","miRNA , mRNA, Methy L  PCA Smote (2000 var revised new grid) "  )) 

# PCA - H methy 2000 vars - 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"3543c2ce-8674-4639-acf8-aa61e518710b -us","miRNA , mRNA, Methy H PCA Undersample (2000 var revised new grid)"  )) 
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"3543c2ce-8674-4639-acf8-aa61e518710b -sm","miRNA , mRNA, Methy H  PCA Smote (2000 var revised new grid) "  )) 

# note 01.12.2020 new uuid added: "038dc8ea-d102-4812-ac65-066e96bd0af6" to list mRNA and Metlylation results X removed...
# 038dc8ea-d102-4812-ac65-066e96bd0af6 (RF)
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"038dc8ea-d102-4812-ac65-066e96bd0af6 -sm","mRNA and Metly  hbybit Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"038dc8ea-d102-4812-ac65-066e96bd0af6 -us","mRNA and Metly  hbybit Smote  sm"  )) 
# 95d3c9a9-a67d-40df-8000-d3428476723e (PCA)
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"95d3c9a9-a67d-40df-8000-d3428476723e -sm","mRNA and Metly  PCA Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"95d3c9a9-a67d-40df-8000-d3428476723e -us","mRNA and Metly  PCA  Smote  sm"  ))

# NOte 02.01.2021- added new uuid:b36fd6ba-c2a2-4b88-8035-0ac21679ac52 for Metlylation results 
# b36fd6ba-c2a2-4b88-8035-0ac21679ac52 (RF)
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"b36fd6ba-c2a2-4b88-8035-0ac21679ac52 -us","Methy only  Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"b36fd6ba-c2a2-4b88-8035-0ac21679ac52 -sm","Methy only  Smote  sm"  )) 
# 036b2a27-1387-4c0c-97f0-245512d953db (PCA) -2000 vars
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"036b2a27-1387-4c0c-97f0-245512d953db -us","Methy only PCA  Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"036b2a27-1387-4c0c-97f0-245512d953db -sm","Methy only PCA Smote  sm"  )) 


# Note :03.01.2021 7f4e41c0-80e6-4932-b16f-34f1dafb8a54 to list mRNA results ( first 500 features)
# 7f4e41c0-80e6-4932-b16f-34f1dafb8a54 (RF)
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"7f4e41c0-80e6-4932-b16f-34f1dafb8a54 -us","mRNA  only  Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"7f4e41c0-80e6-4932-b16f-34f1dafb8a54 -sm","mRNA only  Smote  sm"  ))
# a00c2f50-ab69-40be-af16-b7282c612c38 (PCA)
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"a00c2f50-ab69-40be-af16-b7282c612c38 -us","mRNA  only PCA  Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"a00c2f50-ab69-40be-af16-b7282c612c38 -sm","mRNA only  PCA Smote  sm"  ))

# Note :07.01.2021 ( first 500 features) miRNA & Methy

# b1d0deb3-db25-48de-99de-1437158f159e (RF)
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"b1d0deb3-db25-48de-99de-1437158f159e -us","miRNA-Methy  Hybrit  Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"b1d0deb3-db25-48de-99de-1437158f159e -sm","miRNA-Methy   Hybrit Smote  sm"  ))

# cbc78f99-f3e0-47fa-82f8-e689ab04b260 PCA
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"cbc78f99-f3e0-47fa-82f8-e689ab04b260 -us","miRNA-Methy PCA   Undersample us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"cbc78f99-f3e0-47fa-82f8-e689ab04b260 -sm","miRNA-Methy  PCA Smote  sm"  ))


# 8bd400ed-e910-4f8c-8680-0fddabb3e3f0 3 set retry PCa
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"8bd400ed-e910-4f8c-8680-0fddabb3e3f0 -sm","miRNA-mRNA-Methy PCA   SMOTE us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"8bd400ed-e910-4f8c-8680-0fddabb3e3f0 -us","miRNA-mRNA-Methy PCA Undersample  sm"  ))

# Note :10.01.2021 ( first 500 features) mRNA & Methy 
# 73e8a4e5-0431-47c0-a214-3f5e3b9868d7 - Hybrit
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"73e8a4e5-0431-47c0-a214-3f5e3b9868d7 -sm","mRNA-Methy  Hybrit  Smote us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"73e8a4e5-0431-47c0-a214-3f5e3b9868d7 -us","mRNA-Methy Hybrit Undersample  sm"  ))

# ?? PCA
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"2e156abe-48d4-44ee-8a46-c456939dad3d -sm","mRNA-Methy  PCA  Smote us"  ))
output <- rbind(output, CalculateAccuracyResultsv2(Result_Matrix_April_Rf,"2e156abe-48d4-44ee-8a46-c456939dad3d -us","mRNA-Methy PCA Undersample  sm"  ))


View(output)

logToFile("Predict","training_results_02012021.csv",output,"L")


output <- NULL

# rm(dataSample, dataSamples,dataSampleList, dataSample.1, dataSample.2, dataSample.3, dataSample.4, dataSample.5, dataSample.6, dataSample.7, dataSample.8, dataSample.9, dataSample.10)
# svmLinear,ranger, svmPoly,nnet,naive_bayes, svmRadialGrid,adaBoost

####################BOXPLOTS#################

compositeResults <- rbind(Result_Matrix_April_Rf$result,ResultsPCA$result)
compositeResults[compositeResults$`model_list[j]` %in% "svmRadial",]$`model_list[j]` <- c(rep("svmRadialGrid",30))

trials <- data.frame(model=c("nnet","svmLinear", "svmPoly",  "svmRadialGrid","ranger","adaBoost"),#,"naive_bayes"),
                    names=c("NNET","SVM (Linear)", "SVM(Poly)", "SVM(Radial)","Random Forest","AdaBoost"))#,"Naive Bayes"))

par(mfrow = c(2, 2))

for (i in 1:6) {
    

  df <-
    data.frame(
      "trail_code" = c("0c4aba8a-6506-4ce0-9a94-261934682414 -us",
                       "a00c2f50-ab69-40be-af16-b7282c612c38 -sm",
                       "036b2a27-1387-4c0c-97f0-245512d953db -sm",
                       "ab4cc487-2d45-46eb-9784-4449a731be27 -sm",
                       "73e8a4e5-0431-47c0-a214-3f5e3b9868d7 -sm",
                       "cbc78f99-f3e0-47fa-82f8-e689ab04b260 -sm",
                       "8bd400ed-e910-4f8c-8680-0fddabb3e3f0 -sm",
                       "3ae16ed8-43b1-4d0d-b0ca-08a0ec0f7341 -sm",
                       "3543c2ce-8674-4639-acf8-aa61e518710b -sm"
                      
      ),
      "model_type" = c("ranger","ranger","svmPoly","svmLinear","ranger","svmRadialGrid3","svmLinear","svmLinear","svmLinear"),
     
      "def"= c("1","2","3","4","5","6","7","8","9")
    )
   head(df)
  temp <- CompareAccuracyResult_9Set(compositeResults,df,paste( "Model Comparision For " ,"Biomarker Sets"))
  View(temp)  
}
runcodes <- c("0c4aba8a-6506-4ce0-9a94-261934682414 -us",
              "ab4cc487-2d45-46eb-9784-4449a731be27 -sm",
              "8bd400ed-e910-4f8c-8680-0fddabb3e3f0 -sm",
              "")

par(mfrow = c(1, 1))

names <- c(rep("mRNA",4))

graph_def<-c(", Hbyrid, Undersample (A7)", ", Hbyrid, Smote (B7)", ", PCA, Undersample (C7)",", PCA, Smote (D7)") 

for (i in 1:4) {
  
  df <-
    data.frame(
      "trail_code" = c(
        rep(runcodes[i],6)
      ),
      "model_type" =  trials$model,
      "def"= c(1:6)
    )
  
  temp <- CompareAccuracyResults(compositeResults,df,paste( names[i], graph_def[i]))# compositeResults >
  View(temp)
}

# 6fb58238-0d46-46ff-83e7-b1ef6363e3c9
#####################TUNE##################
par(mfrow = c(2, 4))
df <- Result_Matrix_April_Rf$modelTrainResultsInfo[Result_Matrix_April_Rf$modelTrainResultsInfo$run_uuid %in% "6fb58238-0d46-46ff-83e7-b1ef6363e3c9 -sm",]
dfNnet <- df[df$modelType %in% "nnet",]
dfNnetRun <- dfNnet[dfNnet$run %in% "1",]
dfNnetRun <- dfNnet[dfNnet$run %in% "2",]
dfNnetRun <- dfNnet[dfNnet$run %in% "3",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "4",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "5",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "6",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "7",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "8",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "9",]
tunnePlot(dfNnetRun)
dfNnetRun <- dfNnet[dfNnet$run %in% "10",]
tunnePlot(dfNnetRun)


##################VARIABLE SELECTION#############
# 2 ab4cc487-2d45-46eb-9784-4449a731be27 -sm  svmLinear   2
lvariables <- rbind(pool$miRNAFeatures)
setorder(lvariables,P_Val) 

lvariablesList <- as.vector(lvariables$probName)  
View(lvariablesList)

  
for (runn in c(1:10)) {
  modelInpts <-  fxGetModelDataInputs (dataSample=dataSampleList[[1]],mode="vars",variableList=lvariablesList,runtypes=c("sm"),
                                       usePrev=FALSE,data=NULL,preruuid="")
  
  modelInpts$sm$importance$run_uuid <- c(rep("ab4cc487-2d45-46eb-9784-4449a731be27 -sm", nrow( modelInpts$sm$importance)))
  modelInpts$sm$importance$run <- c(rep(runn, nrow( modelInpts$sm$importance)))
  modelInpts$sm$importance$modelType <- c(rep("common", nrow( modelInpts$sm$importance)))
  colnames(modelInpts$sm$importance) <- colnames(Result_Matrix_April_Rf$importance)
  
  Result_Matrix_April_Rf$importance <- rbind(Result_Matrix_April_Rf$importance, modelInpts$sm$importance)
  
}


importance <- Result_Matrix_April_Rf$importance[Result_Matrix_April_Rf$importance$run_uuid %in% "0c4aba8a-6506-4ce0-9a94-261934682414 -us",]

setorder(importance,importance)

selection <- importance[1:600,]$columnName

selection <- unique( unlist(lapply(selection,  function(x) {str_replace_all(x,"_x_","-")})))


for( i in 1:10)
{
  selection <- calculateCutOf(importance[importance$run %in% i,])
  l <- length(unlist(selection))
  selection <- as.data.frame(selection)
  selection$ds  <- c( rep(i, l ))
  selectedVariables <- rbind(selectedVariables,selection)
}

setorder(importance,importance)

selection <- importance[1:500,]$columnName

#Result_Matrix_April_Rf$importance[Result_Matrix_April_Rf$importance$run_uuid %in% "ab4cc487-2d45-46eb-9784-4449a731be27 -sm",]$run <- c(rep(1, nrow( modelInpts$sm$importance)))

#modelInpts$sm$importance

# 3 1318c3f5-91e5-469b-b230-bb3e8de6e836 -sm  svmLinear   3
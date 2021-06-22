###############################################
# ADI              :  refined Orchestrator 
# YAZAR            : AYSEGUL KUTLAY
# ACIKLAMA         :    
# GUNCELLEME TAR?H?: 06.06.2019
# VERSIYON         : 02
################################################
#load("~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/.RData")

source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/config/environmentVars.R')
source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/config/resources.R')

setwd(global.workingDir)

###### TODO refine data parser as well. // 04.12.2019
######## TODO refine mirtarbase  mappings as well // 04.12.2019 ( output will be used on getCrossTableForSignificancyResults func..)


# call data merge service : miRNA and mRNA expression values will be merged  base on case ic column
dfMerged <- margeData() 
logToFile("Main","MergedData.csv",dfMerged,"D")
# reload
dfMerged <- ReadFromLogFile("Main","MergedData.csv","D")

# Perform normalizarion with center and scale opions. and remove constrant variables.
NormDfMerged <- dataNormalization(3,dfMerged) ## 444 obs 58,980 attr.
NormDfMerged<-removeZeroVariables(NormDfMerged,3) ### 44 obs 58.979 attr.

logToFile("Main","FinalData.csv",NormDfMerged,"D")
# reload ****
NormDfMerged <- read_csv("output/intermediateData/Main_FinalData.csv")

################ TRASFORM AND SPLIT DATA SET ############################

# Split data into classes.****

class1<-NormDfMerged[NormDfMerged$Class %in% "Metastatic", ]  
class2<-NormDfMerged[NormDfMerged$Class %in% "Primary Tumor", ]

class1$class<-rep(c(1),nrow(class1))
class2$class<-rep(c(0),nrow(class2))

# 1: metastastic 0: Primary tumor
class1$class<-factor(class1$class, levels = c(0, 1)) 
class2$class<-factor(class2$class, levels = c(0, 1)) 


######################  Significance test ###################
ttest.results<-significancyTest("ttest",class1[,3:ncol(class1)],class2[3:ncol(class2)])
logToFile("Main","ttest.results.csv",ttest.results,"L")

#reload:
ttest.results <- ReadFromLogFile("Main","ttest.results.csv","L")
ttest.results_detailed <-  getRegulationDetailsResults(Main_ttest_results)
logToFile("Main","ttest.results_detailed.csv",ttest.results_detailed,"L")
View(ttest.results_detailed)
ttest.results_detailed <- ReadFromLogFile("Main","ttest.results_detailed.csv","L")


############# Prepare miRNA referance MAP Data Table
miRTarBase_MTI$hgnc_symbol <- miRTarBase_MTI$`Target Gene`
hgnc_item$hgnc_symbol
logToFile("Main","hgnc_item.csv",hgnc_item,"L")

joind_map <- merge(miRTarBase_MTI,hgnc_item,by="hgnc_symbol")
temp <- joind_map[joind_map$`Species (miRNA)` %in%  "Homo sapiens",]
temp <- TrimMiRNA(temp)

joind_map <- temp
View(joind_map)

logToFile("Main","joind_map.csv",joind_map,"L")
rm(temp)


### Methylation Service need to be run or data neeed to be loaded to run this block

GeneMethylationDetailedResults <- getGeneMethylationDetailedResults(GeneMethylation.ttest.results)
# ############## Commment out on 26 of Feb 2020 /on the fly mapping is beening tried..).################
# rm(class1.train_set,class2.train_set,trainInd,n,class1.test_set,class2.test_set,r,sample,sampleRow,search.results,
#     subset,data.training.pca_smote,data.training_f_rose, data.training_f_smote, data.training_smote,
#     global.type1, global.type2, global.type2_5, global.type2_5_3, global.type2_5_3_4, global.type5, global.type8,
#     global.type8_5_2, global.type_all, genData, genData_2, generatedData,inalDataSet, df, ds, hacide.test, hacide.test,
#     Sonuc_Matrix1, svmLinearGrid, data_example,b,names1, names2)
# 
# rm(data.rose, data.smote, hacide.train, result.part1, result.part2, result.part3, result.part4,result.part5, result.part6, result.part7, result.part8, result.part9, result.part10)
# 
# rm(Sonuc_Matrix, Sonuc_Matrix_Revised, Sonuc_Matrix_try, Sonuc_Matrix_undersample,  Sonuc_Matrix_undersample_p1, Sonuc_Matrix_undersample_2ndlevel,
#    Sonuc_Matrix_undersample_2ndlevel_p1, result1, result2, results, rfrid , selectedVar, selectedVars, svmPolyGrid, svmRadialGrid,
#    Type_1_2, Type_1_3, Type_1_4, Type_1_5, Type_1_6, Type_1_7, Type_1_8, Type_2_3, Type_2_4, Type_2_5,  Type_2_6,  Type_2_7,  Type_2_8)
# 
# rm(trainControl, s,result, Accuracy, AccuracyNull, AccuracyLower, AccuracyPValue, AccuracyUpper, definition.p, definition.p, definition.pcCount, definition.type, Kappa,
#     McnemarPValue, NegPredValue, Prevalence, Sensitivity, Specificity, type1, type2, type5, type8, Type_all,
#    type2,  Type2_5,  Type2_5_3,  Type2_5_3_4,  Type8_5_2)
# # to reload see: logToFile("Main","FinalData.csv",NormDfMerged,"D")
# 
# # SMOTE # commented out 20.12.2019
# 
# # n <- as.data.frame(NormDfMerged)
# #GeneratedData <- SMOTE(n[,-c(1,2,3)],n[,3], K=5)
# 
# # finalDataSet <- as.data.frame( generatedData$data )

###################### 
##################### Perepare DATA SET
finalDataSet <- as.data.frame(NormDfMerged[-c(1,2,(length(NormDfMerged)-2):length(NormDfMerged))]) # remove case id  and class  and undefined columns at  the end 

finalDataSet$class <- NormDfMerged$Class
finalDataSet$case_id <-NormDfMerged$case_id_1

geneMethylationDataSet_norm <- geneMethylationDataSet_norm[-c(1:2)]  # drop class and Case Uuid

finalDataSet <- merge(x=finalDataSet, y=geneMethylationDataSet_norm, by="case_id" )


#finalDataSet <- finalDataSet[ , apply(finalDataSet, 2, function(x) !any(is.na(x)))]  # remove columns with NA

GeneMethylationDetailedResults <- GeneMethylationDetailedResults[GeneMethylationDetailedResults$ProbName %in% colnames(finalDataSet),]
temp <- ttest.results_detailed[ttest.results_detailed$probName %in% colnames(finalDataSet),]
logToFile("Main","finalDataSet.csv",finalDataSet,"D")


###################### Variable SUB Set #########
# Patern.miRNA <-  "NA"
# Patern.mRNA <- "H"
# grpcode <-"6"

###########
##CAN VE RUN AFTER RELEOAD...... 

Grp1.1 <- getPattern("L","Nan","Nan","1.1")
Grp1.2 <- getPattern("H","H","L","2")
Grp2.1 <- getPattern("NA","H","L","3")
Grp2.2 <- getPattern("NA","H","H","3")
Grp3.1 <- getPattern("H","H","L","6")
Grp3.2 <- getPattern("H","H","H","6")

inputSignificanceSet <- rbind(Grp1.1) # ,Grp1.2,Grp2.1,Grp2.2,Grp3.1,Grp3.2)
inputSignificanceSet <- inputSignificanceSet[! is.na( inputSignificanceSet$mRNA_probName), ]
logToFile("Main","inputSignificanceSet.csv",inputSignificanceSet [ order(inputSignificanceSet$miRNA_probName),],"L")


#rm(Grp1,Grp2,Grp3)

selectedVar <- FilterSignificantVariables(inputSignificanceSet,FALSE)

#############################
########Model Train


Sonuc_Matrix_March <-   NULL
modelInfo <- NULL
results <- NULL
p_vals <- c(0.001)
try_count <- 10

pc_count <- c(300)
base_item <- c(0)
model_list<-c("svmLinear","svmRadial" ,"svmPoly","nnet") # rbf eklenecek... "
model_list<-c("svmLinear" ) # rbf eklenecek... "

definition_items <-  c(0) 

definition_list <- c( "All-miRNA:H-L mRNA:H-L Methy: H-L")  # rbf eklenecek...
Grp1 <- getPattern("H","H","H","1.1")

# Grp2 <- getPattern("Nan","H","Nan","2")

inputSignificanceSet <- rbind(Grp1) # ,Grp1.2,Grp2.1,Grp2.2,Grp3.1,Grp3.2)
SelectedVar <- FilterSignificantVariables(inputSignificanceSet,TRUE)

#############################

 Grp1 <-  ttest.results_detailed[ttest.results_detailed$dtype %in% "mRNA",]
 Grp1 <-  Grp1[Grp1$Diff001 %in% "1",]
 Grp1 <-  Grp1[Grp1$direction %in% "L",]

 Grp2 <-  ttest.results_detailed[ttest.results_detailed$dtype %in% "mRNA",]
 Grp2 <-  Grp2[Grp2$Diff001 %in% "1",]
 Grp2 <-  Grp2[Grp2$direction %in% "H",]

 Grp3 <-  ttest.results_detailed[ttest.results_detailed$dtype %in% "miRNA",]
 Grp3 <-  Grp3[Grp3$Diff001 %in% "1",]
 Grp3 <-  Grp3[Grp3$direction %in% "L",]

 Grp4 <-  ttest.results_detailed[ttest.results_detailed$dtype %in% "miRNA",]
 Grp4 <-  Grp4[Grp4$Diff001 %in% "1",]
 Grp4 <-  Grp4[Grp4$direction %in% "H",]

 Grp5 <-  GeneMethylationDetailedResults[GeneMethylationDetailedResults$dtype %in% "Methylation",]
 Grp5 <-  Grp5[Grp5$Diff %in% "1",]
 Grp5 <-  Grp5[Grp5$direction %in% "H",]

 Grp6 <-  GeneMethylationDetailedResults[GeneMethylationDetailedResults$dtype %in% "Methylation",]
 Grp6 <-  Grp6[Grp6$Diff %in% "1",]
 Grp6 <-  Grp6[Grp6$direction %in% "L",]

 colnames(Grp6)<-colnames(Grp1)
 colnames(Grp5)<-colnames(Grp1)
 inputSignificanceSet <- rbind(Grp1,Grp2,Grp3,Grp4) # ,Grp1.2,Grp2.1,Grp2.2,Grp3.1,Grp3.2)

 #SelectedVar <- FilterSignificantVariables(inputSignificanceSet,TRUE)
 SelectedVar <- unique(inputSignificanceSet$probName)
 length(SelectedVar)



##############
run_uuid <- UUIDgenerate(use.time = NA, n=1L)
run_uuid

### Cerate data from selected variables
datasubset <- subset(finalDataSet,select = as.vector( SelectedVar ))
datasubset <- datasubset[ , apply(datasubset, 2, function(x) !any(is.na(x)))]  # remove columns with NA
datasubset$class <- finalDataSet$class

datasubsetclass1<-datasubset[datasubset$class %in% "Metastatic", ]  
datasubsetclass2<-datasubset[datasubset$class %in% "Primary Tumor", ]

datasubsetclass1$class<-rep(c(1),nrow(datasubsetclass1))
datasubsetclass2$class<-rep(c(0),nrow(datasubsetclass2))

# 1: metastastic 0: Primary tumor
datasubsetclass1$class<-factor(datasubsetclass1$class, levels = c(0, 1)) 
datasubsetclass2$class<-factor(datasubsetclass2$class, levels = c(0, 1)) 
datasubset <- rbind(datasubsetclass1,datasubsetclass2)
# 457955ff-9dd3-4c0b-9655-06727372150a All- under sample
# 463d7534-1bed-4452-bb78-aa79f67f8caf mRNA, miRNA- under sample
# e2b454a5-fa50-4ca5-b4ba-c066bb936cc8 All - all
for (i in 1:try_count){
   def <- definition_list[t]
   cat ( paste ( "\n run:", def , "\n"))
    ###############  DATA REORDER & SPLIT ############################
    # class 1: Metastatic
    
   # re order the class 1,create under sample and split from %80 
    class1 <- datasubsetclass1[sample(1:nrow(datasubsetclass1), nrow(datasubsetclass1), replace=FALSE),] 
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
    #data.training <- data.training_f[sample(1:nrow(data.training), nrow(data.training), replace=FALSE),]
    
    #Commet out for Smote
    # generatedData <- SMOTE(data.training[,-c(58978)],data.training[,58978], K=3)
    # data.training <-as.data.frame( generatedData$data )
    
    #bind traning and test set for class 1 and class 2
    data.test <- rbind(class2.test_set,class1.test_set)
    # re order the test set
    #data.test <- data.test[sample(1:nrow(data.test), nrow(data.test), replace=FALSE),]
   
          
    j<-1   
    t<-1
    k<-1
    
   
    
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
    class1_selected <- class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
    model.data.train<-rbind(class1_selected,class2)
    
    class1<-data.test[data.test$class %in% 1, ]  
    class2<-data.test[data.test$class %in% 0, ]
    class1_selected <- class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
    data.test<-rbind(class1_selected,class2)
    test.data_outcomes <- data.test$class
    
    test.data <- predict(data.training.pca, newdata = data.test)
    test.data <- as.data.frame(test.data)
    test.data_outcomes<-data.test$class
    
    XX <- min(pc_count[j], XX)
    
    cat(paste("\n Component Count:", XX,"\n"))
    
    Sonuc_Matrix_try<- model_comp_pca( model.data.train[,1:XX], test.data[,1:XX],test.data_outcomes , model_list)
    Sonuc_Matrix_try$results
    cat("\n")
    print( ( Sonuc_Matrix_try$results[,c(1,2,12,17)])) 
    Sonuc_Matrix_try$results$definition.type  <- rep(c(definition_list[t]), length(length(Sonuc_Matrix_try$results[1])))
    Sonuc_Matrix_try$results$definition.p <-  rep(p_vals[k], length(length(Sonuc_Matrix_try$results[1])))
    Sonuc_Matrix_try$results$definition.pcCount <- rep(XX, length(length(Sonuc_Matrix_try$results[1])))
    Sonuc_Matrix_try$results$trialCode <- rep(run_uuid, length(length(Sonuc_Matrix_try$results[1]))) 
    
    
    results <-  rbind (results,Sonuc_Matrix_try$results)
    modelInfo<- rbind (modelInfo,Sonuc_Matrix_try$model)
    
    Sonuc_Matrix_March$modelInfo <- modelInfo
    Sonuc_Matrix_March$results <- results 
                 
              

}

output <- CalculateAccuracyResults(Sonuc_Matrix_March,"463d7534-1bed-4452-bb78-aa79f67f8caf","miRNA, mRNA, Methylation- undersample"  )
output <- rbind(output, CalculateAccuracyResults(Sonuc_Matrix_March,"b44b6a2f-7af7-4d3a-9c21-0402297ae7f2","miRNA, mRNA- undersample"  )) 
output <- rbind(output, CalculateAccuracyResults(Sonuc_Matrix_March,"e2b454a5-fa50-4ca5-b4ba-c066bb936cc8","miRNA, mRNA, Methylation- All Samples"  )) 
output <- rbind(output, CalculateAccuracyResults(Sonuc_Matrix_March,"b04a0120-8bb8-4543-86a3-0f6c1b5dc0c7","miRNA, mRNA- All Samples"  )) 


View(output)

#Sonuc_Matrix_March$results$definition.type[603:642] <- rep("miRNA:NA mRNA:L", 40)

View(Sonuc_Matrix_March$results) 
logToFile("RunModel",paste("Sonuc_Matrix_March.csv"), Sonuc_Matrix_March ,"L")

save(Sonuc_Matrix_March,file="Sonuc_Matrix_March.RData")
march2 <- load(file="~/Sonuc_Matrix_March.RData")

Sonuc_Matrix_March$modelInfo <- rbind(modelInfo,Sonuc_Matrix_March$modelInfo )
Sonuc_Matrix_March$results <- rbind(results ,Sonuc_Matrix_March$results)

Sonuc_Matrix_March_backup <-Sonuc_Matrix_March
save(Sonuc_Matrix_March_backup,file="Sonuc_Matrix_March_backup.RData")

rm(class1, class1.test_set, class1_selected, class1.train_set, class2, class2.test_set, class2.train_set, data.test, data.test.pca, data.test_f,
   data.training, data.training.pca, data.training_f,Sonuc_Matrix_try, model.data.train)







# for rest see analysis folder...

############  LOG
# 
# for (i in 1:try_count){
#   ###############  DATA REORDER & SPLIT ############################
#   # class 1: Metastatic
#   # re order the class 1 and split from %70
#   class1 <- class1[sample(1:nrow(class1), nrow(class1), replace=FALSE),]
#   # class1 <-  class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
#   trainInd <- round(nrow(class1)*0.8) 
#   # create test and train sets for class 1
#   class1.train_set <- class1[1:trainInd,]
#   class1.test_set <- class1[-(1:trainInd),]
#   
#   #reorder the class 2 and split from %80
#   class2 <- class2[sample(1:nrow(class2), nrow(class2), replace=FALSE),]
#   
#   trainInd <- round(nrow(class2)*0.8) 
#   # create test and train sets for class 2
#   class2.train_set<-class2[1:trainInd,]
#   class2.test_set <- class2[-(1:trainInd),]
#   
#   #bind traning sets for class 1 and class 2
#   data.training_f<-rbind(class1.train_set,class2.train_set)
#   #re order the sets
#   data.training_f <- data.training_f[sample(1:nrow(data.training_f), nrow(data.training_f), replace=FALSE),]
#   
#   #Commet out for Smote
#   # generatedData <- SMOTE(data.training_f[,-c(58978)],data.training_f[,58978], K=3)
#   # data.training_f <-as.data.frame( generatedData$data )
#   
#   #bind traning and test set for class 1 and class 2
#   data.test_f<- rbind(class2.test_set,class1.test_set)
#   # re order the test set
#   data.test_f <- data.test_f[sample(1:nrow(data.test_f), nrow(data.test_f), replace=FALSE),]
#   
#   for (t in 1:length( definition_list)) {
#     
#     for (k in 1:length( p_vals)) {
#       
#       
#       typeList <- c( base_item,  definition_items[t])
#       cat( "\n items: ")
#       print(typeList)
#       ### selectedVar <- getSignificantVariables(flaggedSigResult, p_vals[k],typeList) ##### the logic have been modfid ne method is written. on 26th of FEB
#       
#       
#       # Prepare Sub Set of significant variables....
#       data.training <- subset(data.training_f,select = as.vector( SelectedVar ))
#       data.training$class <- data.training_f$class
#       
#       data.test <- subset(data.test_f,select = as.vector( SelectedVar ))
#       data.test$class <- data.test_f$class
#       
#       
#       for (j in 1:length( pc_count)) {
#         def <- paste(definition_list[t],p_vals[k],pc_count[j], sep = "_")
#         cat ( paste ( "\n run:", def , "\n"))
#         
#         ################ PCA ############################
#         # principal component Analaysis 
#         
#         data.training.pca<- prcomp(data.training[,-c(ncol(data.training))])
#         
#         # Comment out for PCA Logs and figures..
#         #compute standard deviation of each principal component
#         std_dev <- data.training.pca$sdev
#         #compute variance
#         pr_var <- std_dev^2
#         #check variance of first 100 components
#         # pr_var[1:100]
#         #proportion of variance explained
#         prop_varex <- pr_var/sum(pr_var)
#         sumcomp <-0
#         for (comp in 1: length( prop_varex)) {
#           sumcomp <- sumcomp + prop_varex[comp]
#           if (sumcomp < 0.98) { XX <- comp }
#         }
#         #scree plots
#         # filePart=paste(global.workingDir,"/output/figures/","PCA","_", sep = "")
#         # png(filename=paste(filePart,"ProportionVariancExplaine_", def,"png"))
#         # plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained",  type = "b")
#         # dev.off()
#         # png(filename=paste(filePart,"CumProportionVariancExplaine_", def,"png"))
#         # plot(cumsum(prop_varex), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained",  type = "b")
#         # dev.off()
#         
#         # Prepare for prediction
#         data.test.pca <- predict(data.training.pca, newdata = data.test)
#         
#         
#         ############# MODEL RUN ##############
#         
#         cat(paste("Try:",i, "Set:", t, "Metod:", k, j,  "\n"))
#         model.data.train <- data.frame(Class = data.training$class , data.training.pca$x)
#         #transform test into PCA
#         
#         # re order dataset
#         model.data.train  <- model.data.train[sample(1:nrow(model.data.train), nrow(model.data.train), replace=FALSE),]
#         data.test  <- data.test[sample(1:nrow(data.test), nrow(data.test), replace=FALSE),]
#         
#         # under sampling for metastastic to generate new set 
#         class1<-model.data.train[model.data.train$Class %in% 1, ]  
#         class2<-model.data.train[model.data.train$Class %in% 0, ]
#         class1_selected <- class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
#         model.data.train<-rbind(class1_selected,class2)
#         
#         class1<-data.test[data.test$class %in% 1, ]  
#         class2<-data.test[data.test$class %in% 0, ]
#         class1_selected <- class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
#         data.test<-rbind(class1_selected,class2)
#         test.data_outcomes <- data.test$class
#         
#         test.data <- predict(data.training.pca, newdata = data.test)
#         test.data <- as.data.frame(test.data)
#         test.data_outcomes<-data.test$class
#         
#         XX <- min(pc_count[j], XX)
#         
#         cat(paste("\n Component Count:", XX,"\n"))
#         
#         Sonuc_Matrix_try<- model_comp_pca( model.data.train[,1:XX], test.data[,1:XX],test.data_outcomes , model_list)
#         Sonuc_Matrix_try$results
#         cat("\n")
#         print( ( Sonuc_Matrix_try$results[,c(1,2,12,17)])) 
#         Sonuc_Matrix_try$results$definition.type  <- rep(c(definition_list[t]), length(length(Sonuc_Matrix_try$results[1])))
#         Sonuc_Matrix_try$results$definition.p <-  rep(p_vals[k], length(length(Sonuc_Matrix_try$results[1])))
#         Sonuc_Matrix_try$results$definition.pcCount <- rep(XX, length(length(Sonuc_Matrix_try$results[1])))
#         Sonuc_Matrix_try$results$trialCode <- rep(run_uuid, length(length(Sonuc_Matrix_try$results[1]))) 
#         
#         
#         results <-  rbind (results,Sonuc_Matrix_try$results)
#         modelInfo<- rbind (modelInfo,Sonuc_Matrix_try$model)
#         
#         Sonuc_Matrix_March$modelInfo <- modelInfo
#         Sonuc_Matrix_March$results <- results 
#         
#       }# end for ppc count for
#     } # end for definition_list
#   } # end  for p_vals
#   
# }




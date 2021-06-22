###############################################
# ADI              :  Data Prepare Util Functions script
# YAZAR            : AYSEGUL KUTLAY
# ACIKLAMA         :    
# GUNCELLEME TAR?H?: 06.06.2019
# VERSIYON         : 02
################################################
# performs normalization  operation for gievn data accoridng to given model
doNormalization <-function(Data,Model)
{
  
  min_var=min(Data)
  max_var=max(Data)
  sdv_var=sd(Data)
  mean_var=mean(Data)
  
  normD<-NULL
  
  l=length(Data)
  
  if(Model=="minmax")
  {
    MINMAXNORM_VAR126=Data# min-max26
    for (i in 1:l)  MINMAXNORM_VAR126[i] <- ((Data[i]-(min_var-0.000001))/(max_var-(min_var-0.000001)))
    normD<-MINMAXNORM_VAR126
  }
  
  if(Model=="zscore")
  {
    ZSCORENORM_VAR126=Data#z score
    for (i in 1:l)  ZSCORENORM_VAR126[i] <- ((Data[i]-mean_var)/(sdv_var))
    normD<-ZSCORENORM_VAR126
  }
  
  if(Model=="log")
  {
    
    normD<-as.numeric(unlist(log(Data+transform(rep(c(0.0000001),NROW(Data))))))
  }
  
  
  normD
}
prepareData <-function()
{
  data<-read.csv2("1495397748032-manifest.csv")
  myMetaData<-as.data.frame(data)
  
  for (j in 1:nrow(myMetaData)-1) # ncol(myMetaData)-1
  {
    filename<- as.character( myMetaData$name[j])
    
    
    if(length(grep( "mirnas.quantification", filename)) >0){
      dataFile<-read.delim(filename)
      miRNA_Data<-t(as.vector( dataFile$reads_per_million_miRNA_mapped))
      colnames(miRNA_Data)<-as.character( dataFile$miRNA_ID)
      d<-cbind(myMetaData[j,],miRNA_Data)
      #as.vector( myMetaData[2,])
      #as.vector(miRNA_Data)
      if(j>2)
      {
        miRNA_Data_Segment<-rbind(miRNA_Data_Segment,d)
      }else{
        miRNA_Data_Segment<-d
      }
    }
    
  }#end for
  
  write.csv(file="melonomaData.csv",miRNA_Data_Segment)
  miRNA_Data_Segment
}

performChi2Test <-function(myData)
{
  results.name<-c()
  results.nameCrossvar<-c()
  results.significant<-c()
  t<-1
  
  for (j in 4:ncol(myData)-1)#ncol(data))
  {
    
    if(sum(is.na(myData[,j]))>100)
    {
      d<-cbind(myData[,j],myData[,ncol(myData)])
      chi2.result<-chisq.test(na.omit(d))
      # t.result= t.test(myData.Miss[,126],data.label, var.equal=FALSE, paired=FALSE)
      # chi2.result#$$p-value
      if(chi2.result$p.value>0.05)     
        results.significant[t] =1  
      else
        results.significant[t] =0
      
      results.name[t]<-names(myData[j])
    }
    #  results.nameCrossvar[t]<-names(data[j])
    t=t+1
  }
  
  results.all<-data.frame(results.name,results.significant)
  write.csv(file="significanceResultsWithChi2.csv",results.all)
  results.all
}


significancyTest<-function(type,class1,class2){
  size <- ncol(class1)
  cat(size)
  result <- data.frame(ProbName= rep(NA,size),Ftest=rep(NA,size),Ttest=rep(NA,size),Diff=rep(NA,size),P_Val=rep(NA,size), P_Val_Less=rep(NA,size),P_Val_Greater=rep(NA,size))
 
  
   for (i in 3:(ncol(class1)-1))
   {
     
    cat(paste("\n", i) )
    
    if (type=="ttest")
    {
        propResult<-TTest(colnames(class1[i]), class1[i],  class2[i])
    }else if(type=='pearson'){
      propResult<-pearsonTest(colnames(class1[i]), class1[i],  class2[i])
    
      
    } 
    
    result$ProbName[i] <-  as.character( propResult$probName[1])
     result$Ftest[i] <- as.character( propResult$Ftest[1])
     result$Ttest[i]  <- as.character(propResult$Ttest[1])
     result$Diff[i]  <-  propResult$Diff[1]
     result$P_Val[i]  <- propResult$P_Val[1]
     result$P_Val_Less[i]  <- propResult$P_Val_Less[1]
     result$P_Val_Greater[i]  <- propResult$P_Val_Greater[1]
     
      
     cat("\014") 
   }
  
  result
  
}


library(MASS) 
library(stats)
# read csv files
pearsonTest <-function(parameterName,DataX,DataY)
{
  
  alphaTest=0.99
  
  prob=parameterName
  
  Ftest=c()
  Ttest=c()
  Diff=c()
  P_Val=NULL
  i=1
  x<-DataX[,1]# as.vector(sapply(DataX, as.numeric))
  y<-DataY[,1]#as.vector(sapply(DataY, as.numeric))
  fVar<-cor.test(x,y,alternative="two.sided",method="pearson")
  fVar$p.value
  
  P_Val=fVar$p.value
  if(fVar$p.value>alphaTest)
  {
    Ttest=(' differentially expresse between Two groups')     
    Diff=1
    
  }else
  {
    Ttest=('Not Differentially expresse between Two groups')     
    Diff=0
  }
  
  
  
  GeneTestData<-data.frame(prob,Ttest,Diff,P_Val)
  View(GeneTestData)
  GeneTestData
}

library(MASS) 
library(stats)
# read csv files
TTest <-function(probName,DataX,DataY)
{
  
  # DataX<-class1[i]
  # DataY<-class2[i]
 # DataX
  # DataY
  #threasholds
  test.alphaVar=0.001
  test.alphap=0.05
  test.conflevel=(1-test.alphap)
  test.variance=FALSE
  
  
  #result Variables
  Ftest=c()
  Ttest=c()
  Diff=c()
  P_Val=NULL
  P_Val_Less=1
  P_Val_Greater=1
  test.variance=FALSE
  
  
  #Start test
  Grp1<- DataX[,1];
  Grp2<-  DataY[,1];
  fVar=var.test(Grp1, Grp2) 
  fVar
  if(fVar$p.value>test.alphaVar)
  {
    Ftest='Variances are homogeneous'
    test.variance=TRUE
    
  }else
  {
    Ftest=('Variances are Not homogeneous')
    test.variance=FALSE
    
  }
  
  t.result=t.test(Grp1, Grp2 , alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = test.variance, conf.level = test.conflevel)
  
  
  P_Val=t.result$p.value
  if (test.alphap<P_Val) 
  {
    Ttest=('Not differentially expressed between Two groups')     
    Diff=0
  }else{
    
    Ttest=('Differentially expressed between Two groups')     
    Diff=1 
    
     t.less=t.test(Grp1, Grp2 , alternative = c( "less"), mu = 0, paired = FALSE, var.equal =  test.variance, conf.level = test.conflevel)
     t.greater=t.test(Grp1, Grp2 , alternative = c("greater"), mu = 0, paired = FALSE, var.equal = test.variance, conf.level = test.conflevel)
     P_Val_Less<-t.less$p.value
     P_Val_Greater<-t.greater$p.value
  }
  
  
  GeneTestData<-data.frame(probName,Ftest,Ttest,Diff,P_Val, P_Val_Less,P_Val_Greater)
   
  GeneTestData
}


variableSelection<-function(intputData,filename){
  # intputData <- htseq_data_step_one
  # filename<-"mrna"
  intputData<-intputData[intputData$sample_type %in% c("Metastatic","Primary Tumor"),]
  #remore metadata
  
  dataSegment<-intputData[c(24:ncol(intputData),9)]
  cat(paste("/n Prepare For Test",ncol(dataSegment)))
  for (i in 1:(ncol(dataSegment)-1)){
    cat(paste("\n", i))
    dataSegment[i]<-as.double( unlist(dataSegment[i]) )
    cat("\014")  
  }
  
  #omit vars with zero entry
  remove<-c()
  cat(paste("/n  omit null values Total:",ncol(dataSegment)))
  for(i in 1:(ncol(dataSegment)-1))
  {
    cat(paste("\n", i))
    if(min(dataSegment[i])== max(dataSegment[i])){
      remove<-cbind(remove,(i))
    }
    cat("\014")  
  }
  dataSegment2<-dataSegment
  #subract
  if (length(remove)>0)
  {
    a<-dataSegment[, as.vector( remove)]
    dataSegment<-dataSegment[, -as.vector( remove)]
  }
  
  
  # significancy test
  dataSegment$sample_type<-as.factor(dataSegment$sample_type)
  finnalData<-dataSegment
  class1<-finnalData[finnalData$sample_type %in% "Metastatic", ]  
  class2<-finnalData[finnalData$sample_type %in% "Primary Tumor", ]
  
  #Chi2Test<-performChi2Test(finnalData)
  ttestResult <-significancyTest('ttest',class1,class2)   
  
  #pearsonResults<-significancyTest('pearson',class1,class2)
  
  logToFile("libs_variableSelection",paste(filename,"_significancyTestResultst.csv",sep=""),ttestResult,"L")
  
  # write.csv(ttestResult,paste(filename,"_significancyTestResultst.csv",sep=""))
  
  significantGeneSet<-ttestResult[ttestResult$Diff %in% "1",]
  
  selected<-c()
  cat(paste("/n  Start Testing Total:",ncol(finnalData)))
  for(i in 1:(ncol(finnalData)-1))
  { 
    cat(paste("\n", i))
    geneName<-colnames(finnalData[i])
    checkResult<-significantGeneSet[significantGeneSet$probName %in% geneName,]
    if(nrow(checkResult)>0){
      selected<-cbind(selected,(i))
    }
    cat("\014")  
  }
  selected<-cbind(selected,ncol(finnalData))
  
  finalData_selected_geneSet<-finnalData[, as.vector( selected)]
  
  logToFile("libs_variableSelection",paste(filename,"_finalData_selected_geneSet.csv",sep=""),finalData_selected_geneSet,"L")
  
  # write.csv(finalData_selected_geneSet,paste(filename,"_finalData_selected_geneSet.csv",sep=""))
  
  finalData_selected_geneSet
  
  
}

runModel<-function(class1,class2,methodName){
  
  # methodName<-"rf"
  # class1_selected <- class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
  data<-rbind(class1,class2)
  data$Class<-factor(data$Class, levels = c(0, 1)) 
  data <- data[sample(1:nrow(data), nrow(data), replace=FALSE),]
  n<-nrow(data)
  
  Model_Data<-data[sample(nrow(data),n),] 
  trainInd <- round(nrow(Model_Data)*0.7) 
  train_set<-Model_Data[1:trainInd,]
  test_Model_Data_outcome <- Model_Data[-(1:trainInd),1]
  #   test_Model_Data_outcome<-factor(test_Model_Data_outcome$Class, levels = c(0, 1)) 
  test_set <-as.data.frame( Model_Data[-(1:trainInd),-1])
  
  model<-train(Class~.,
               data=train_set, 
               method=methodName,
               preProcess = c("center","scale"),
               #  tuneGrid=svmGrid,
               trControl = trainControl(method = "cv",number=10))
  
  
  model
}

model_comp<-function(class1,class2,model_list,times){
  
  results<-data.frame(Sensitivity<-numeric(0),Specificity<-numeric(0) ,  PosPredValue<-numeric(0) ,
                      NegPredValue<-numeric(0) ,      Prevalence<-numeric(0) , DetectionRate<-numeric(0) , 
                      DetectionPrevalence <-numeric(0) ,BalancedAccuracy<-numeric(0) ,
                      Accuracy  <-numeric(0) ,   Kappa <-numeric(0) ,AccuracyLower<-numeric(0) ,
                      AccuracyUpper<-numeric(0) , AccuracyNull<-numeric(0) , AccuracyPValue<-numeric(0) , McnemarPValue<-numeric(0),model<-character(0))
  
  for (i in 1:times)
  {
    class1_selected <- class1[sample(1:nrow(class1), nrow(class2), replace=FALSE),]
    data<-rbind(class1_selected,class2)
    data$Class<-factor(data$Class, levels = c(0, 1)) 
    #satirla?? kar????t??r30
    data <- data[sample(1:nrow(data), nrow(data), replace=FALSE),]
    n<-nrow(data)
    
    Model_Data<-data[sample(nrow(data),n),] 
    trainInd <- round(nrow(Model_Data)*0.7) 
    train_set<-Model_Data[1:trainInd,]
    test_Model_Data_outcome <- Model_Data[-(1:trainInd),1]
    test_Model_Data_outcome<-factor(test_Model_Data_outcome, levels = c(0, 1)) 
    
    test_set <-as.data.frame( Model_Data[-(1:trainInd),-1])
    cat(i)
    
    results<- trainModels(train_set,test_set,model_list,test_Model_Data_outcome) 
    
  }
  results
}

model_comp_pca<-function(train_set,test_set,test_Model_Data_outcome,model_list){
  
  
  results<- trainModels(train_set,test_set,model_list,test_Model_Data_outcome) 
  results
}

trainModels<-function(train_set,test_set,model_list,test_Model_Data_outcome){
  
    results<-data.frame(Sensitivity<-numeric(0),Specificity<-numeric(0) ,  PosPredValue<-numeric(0) ,
                      NegPredValue<-numeric(0) ,      Prevalence<-numeric(0) , DetectionRate<-numeric(0) , 
                      DetectionPrevalence <-numeric(0) ,BalancedAccuracy<-numeric(0) ,
                      Accuracy  <-numeric(0) ,   Kappa <-numeric(0) ,AccuracyLower<-numeric(0) ,
                      AccuracyUpper<-numeric(0) , AccuracyNull<-numeric(0) , AccuracyPValue<-numeric(0) , 
                      McnemarPValue<-numeric(0),model<-character(0), definition.type <-character(0), 
                      definition.p <- character(0) , definition.pcCount <- character(0))
  
  response <- NULL
  models <-NULL
  model <- NULL
  bestTunes <- NULL
  
  for (j in 1:length(model_list))
  {
    
   # #tunnig parameters....
    # svmRadialGrid <- expand.grid(sigma= 2^c(-25:15), C= 10^(-4) * (10:50))  >>  C= 10^(-2) * (2:20)) 
    #svmRadialGrid <- expand.grid(sigma= 2^c(-15:5),C= 10^(-4) * (20:50) )
    # US
    svmRadialGrid <- expand.grid(sigma= 2^c(-15:-5),C= 10^(-2) * (10:25) )  
    svmRadialGrid1 <- expand.grid(sigma= 2^c(-7),C= 10^(-4) * (1) )  
    # SM
    svmRadialGrid2 <- expand.grid(sigma= 2^c(-15:-5),C=  5 * c( 10^(-4) ) )  # 2-20
    svmRadialGrid3 <- expand.grid(sigma= 2^c(-15:-5),C=  10^(-4) * c(5,10,15,20,25) )  # 2-20
    
    svmLinearGrid <- expand.grid( C= 10^(-4) * (20:150))
    svmPolyGrid<- expand.grid(.degree=(2:4), .scale= 0.01, .C= 2^c(-10:-5) )# 2^c(-5:1))  # c:-10,1, c=(-5,5)
    svmPolyGrid2<- expand.grid(.degree=(6), .scale= 0.001, .C= 2^c(-10) )# 2^c(-5:1))  # c:-10,1, c=(-5,5)
    
    dnngrid <- expand.grid(layer1 = 1:25, layer2 = 0:25, layer3 = 0:25)
    dnngrid$hidden_dropout <- 0
    dnngrid$visible_dropout <- 0
    
    nngrid=expand.grid(.size=(1:3),.decay=c(0.01,0.1,1,3,5,10,20))
    
    nngrid2=expand.grid(.size=(1:8),.decay=c(0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.04,0.06,0.08,0.1,1,3,5,10))
   
     nngridrf3selected=expand.grid(.size=(2),.decay=c(1))
    
    
     rfgrid <- expand.grid(
      .mtry = 4:6,
      .splitrule = c("gini"),
      .min.node.size = c(4:10 )
    )
    
    adaboostgrid <- expand.grid(mfinal = (4:9)*2, maxdepth = c(3, 10), coeflearn = c("Breiman")) # c("Breiman", "Freund", "Zhu")
    
    adaboostgrid2 <- expand.grid(mfinal = 2^(1:4), maxdepth = c(2, 6), coeflearn = c("Breiman")) # c("Breiman", "Freund", "Zhu")
    
    maxnetGrid <- expand.grid(
      .layer1 = (1:6),
      .layer2 = (1:6),
      .layer3 = (1:6),
      .learning.rate=10^(-2) * (10:50),
      .momentum=10^(-3),
      .dropout=10^(-3),
      .activation="sigmoid"
    )
    
    nbgrid <- expand.grid(
      .laplace= 10^(-2)*(2:10),
      .usekernel=TRUE,
      .adjust=10^(-1) * c(5,10,15,20)
    )
    
    nbgridsm <- expand.grid(
      .laplace= 10^(-2)*c(2,4,6,8,10),
      .usekernel=TRUE,
      .adjust=10^(-2) * c(5,10,15,20)
    )
    
    nbgridus2 <- expand.grid(
      .laplace= 10^(-3) * c(10:15),
      .usekernel=TRUE,
      .adjust=10^(-2) * c(10:15)
    )
    
    trainControl <-trainControl(method = "repeatedcv", number = 10, repeats = 5) 
    
    cat(paste("\n ",model_list[j]))
    
    #Modelleri kalistiriyoruz. Normalizasyon ve 10 fold cross valdation uyguluyoruz  see train control...
    if("svmRadialGrid"== model_list[j]) 
    {
      cat(" Started ")
      model<-train(class~.,
                   data=train_set, 
                   method="svmRadial",
                   preProcess = c("center","scale"),
                   tuneGrid=svmRadialGrid, #*********
                   trControl = trainControl)
    } else if("svmRadialGrid3"== model_list[j]) 
    {
      cat(" Started.. ")
      model<-train(class~.,
                   data=train_set, 
                   method="svmRadial",
                   preProcess = c("center","scale"),
                   tuneGrid=svmRadialGrid3, #********
                   trControl = trainControl)
    } else if("svmLinear"== model_list[j]) {
      model<- train(class~.,
                   data=train_set, 
                   method=model_list[j],
                   preProcess = c("center","scale"),
                   tuneGrid=svmLinearGrid,
                   trControl = trainControl)
      
    }else if("svmPoly"== model_list[j]){
    
      model<-train(class~.,
                   data=train_set, 
                   method="svmPoly",
                   preProcess = c("center","scale"),
                   tuneGrid=svmPolyGrid2,
                   trControl = trainControl)
    }else if("dnn"== model_list[j]){
      model<-train(class~.,
                   data=train_set, 
                   method=model_list[j],
                   preProcess = c("center","scale"),
                   tuneGrid=dnngrid,
                   trControl = trainControl)
    }else if("nnet"== model_list[j]){
      model<-train(class~.,
                   data=train_set, 
                   method="nnet",
                   preProcess = c("center","scale"),
                   tuneGrid=nngridrf3selected,
                   trControl = trainControl)
    }else if("xyf"== model_list[j]){
      model<-train(class~.,
                   data=train_set, 
                   method=model_list[j],
                   preProcess = c("center","scale"),
                   trControl = trainControl)
    }else if("ranger"== model_list[j]){
      model<-train(class~.,
                   data=train_set, 
                   method=model_list[j],
                   tuneGrid=rfgrid,
                   preProcess = c("center","scale"),
                   trControl = trainControl)
    }else if("naive_bayes2"== model_list[j]){
      model<-train(class~.,
                   data=train_set, 
                   method="naive_bayes",
                   tuneGrid=nbgridus2 , #*******
                   preProcess = c("center","scale"),
                   trControl = trainControl)
    }else if("adaBoost"== model_list[j]){
             model<-train(class~.,
                    data=train_set, 
                    method = "AdaBoost.M1", 
                    trControl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = FALSE) ,
                    tuneGrid = adaboostgrid,
                    preProcess = c("center", "scale"))
    }else if("mxnet"== model_list[j]){
      model<-train(class~.,
                   data=train_set, 
                   method = "mxnet", 
                   trControl = trainControl,
                   tuneGrid = maxnetGrid,
                   preProcess = c("center", "scale"))
    }
    else{
      model<-train(class~.,
                   data=train_set, 
                   method=model_list[j],
                   preProcess = c("center","scale"),
                   tuneLength=150,
                   trControl = trainControl)
      
      
      
    }
    
    prediction <- predict(model, newdata=test_set)
    con_mat<-confusionMatrix(data=prediction, test_Model_Data_outcome)
    
    #con_mat icerisinden gerekli sonuc parametreleri alinir.
    result1<-transpose(as.data.frame(con_mat$byClass))
    colnames(result1)<-rownames(as.data.frame(con_mat$byClass))
    
    result2<-transpose(as.data.frame(con_mat$overall))
    colnames(result2)<-rownames(as.data.frame(con_mat$overall))
    
    result<-cbind(result1,result2,model_list[j])
    
    bestTune <- cbind(model$bestTune,max(na.omit( model$results$Accuracy)))
    bestTunes <-rbind(bestTunes,bestTune)
    results<-rbind(results,result)
    models <- rbind(models, model$results)
  }
  
  
  
  response$results <- results
  response$prediction <-prediction
  response$model <- models
  response$bestTunes <- bestTunes
  response
}

####################################################
#Sensitivity, specificity, accuracy degerlerinin ortalamasini hesaplayan fonksiyon
####################################################
model_perf_comp<-function(a,model_list) {
  
  r<-list()
  for (i in 1:length(model_list))
  {
    res<-subset(a,a[,ncol(a)]==model_list[i])
    r<-list.append(r,res)
  }
  
  results_mean<-data.frame(Model<-character(0),Sensitivity<-numeric(0),Specificity<-numeric(0),Accuracy<-numeric(0),AccuracyPValue<-numeric(0),Kappa<-numeric(0))
  
  for (i in 1:length(r))
  {
    sen_mean<-mean(r[[i]]$Sensitivity)
    sp_mean<-mean(r[[i]]$Specificity)
    ac_mean<-mean(r[[i]]$Accuracy)
    p_Val<-mean(r[[i]]$AccuracyPValue)
    kappa<-mean(r[[i]]$Kappa)
    res_mean<-cbind(model_list[i],sen_mean,sp_mean,ac_mean,p_Val,kappa)
    results_mean<-rbind(results_mean,res_mean)
  }
  results_mean
}

#####################################################
#Modeller sonu?larinda aldigimiz predictionlari dataframe'e atan fonksiyon
#########################################################

model_comp_pred<-function(data,model_list,n){
  
  Model_Data<-data[sample(nrow(data),n),] 
  trainInd <- round(nrow(Model_Data)*0.7) 
  tr_Model_Data<-Model_Data[1:trainInd,]
  test_Model_Data_outcome <- Model_Data[-(1:trainInd),ncol(data)]
  test_Model_Data_Tani <- Model_Data[-(1:trainInd),-ncol(data)]
  pred<-as.data.frame(matrix(nrow=nrow(test_Model_Data_Tani),ncol=length(model_list)))
  
  for (j in 1:length(model_list))
  {
    #Modelleri ?alistiriyoruz. Normalizasyon ve 10 fold cross valdation uyguluyoruz
    
    model<-train(Class~.,data=tr_Model_Data, method=model_list[j],
                 #preProcess = c("center","scale"),
                 trControl = trainControl(method = "cv",number=10))
    prediction <- predict(model, newdata=test_Model_Data_Tani)
    #con_mat<-confusionMatrix(data=prediction, test_Model_Data_outcome)
    pred[,j]<-prediction
  }
  pred<-cbind(pred,test_Model_Data_outcome)
  col_name<-append(model_list,"Class")
  colnames(pred)<-col_name
  pred$Class<-as.factor(pred$Class)
  pred
}


removeZeroVariables <- function(dataSegment,startindex){
  remove<-c()
  i<-startindex
  for (i in startindex:ncol(dataSegment)){
    cat(paste("\n", i))
    dataSegment[i]< as.integer(unlist(dataSegment[i]))
    if(  min(na.omit( dataSegment[i]))==max(na.omit(dataSegment[i]))){
      remove<-cbind(remove,(i))
      
    }
    cat("\014")
  }
  
  
  #subract
  
  dataSegment_2<-dataSegment[ ,-as.vector( remove)]
  dataSegment_2
}


removeNullClass <- function(dataSegment){
  #dataSegment<-miRNA_data_final
  remove<-c()
  for (i in 1:nrow(dataSegment)){
    if( is.na(dataSegment$sample_type[i])){
      remove<-cbind(remove,(i))
    }
  }
  
  #subract
  a<-dataSegment[ as.vector( remove),]
  
  dataSegment_2<-dataSegment[ -as.vector( remove),]
  dataSegment_2
}

CalculateAccuracyResults <- function(accuracy_result, trail_code,GrfgDefinition)
{
  
  accuracy_result <- accuracy_result$results[accuracy_result$results$trialCode %in% trail_code ,]
  nnetResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "nnet", ]  
  
  svmLResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "svmLinear", ]  
  svmPResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "svmPoly", ]  
  svmRResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "svmRadial", ]  
  rfResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "rf", ] 
  
  output <- as.data.frame(  cbind(paste(GrfgDefinition, " nnet"), CI(nnetResults$Sensitivity, ci=0.95)[2], CI(nnetResults$Specificity, ci=0.95)[2], CI(nnetResults$AccuracyPValue, ci=0.95)[2],CI(nnetResults$Accuracy, ci=0.95)[2]))
  output <- rbind(output, as.data.frame(  cbind(paste(GrfgDefinition, " Linear"), CI(svmLResults$Sensitivity, ci=0.95)[2], CI(svmLResults$Specificity, ci=0.95)[2], CI(svmLResults$AccuracyPValue, ci=0.95)[2],CI(svmLResults$Accuracy, ci=0.95)[2])))
  output <- rbind(output, as.data.frame(  cbind(paste(GrfgDefinition, " Poly"), CI(svmPResults$Sensitivity, ci=0.95)[2], CI(svmPResults$Specificity, ci=0.95)[2], CI(svmPResults$AccuracyPValue, ci=0.95)[2],CI(svmPResults$Accuracy, ci=0.95)[2])))
  output <- rbind(output, as.data.frame(  cbind(paste(GrfgDefinition, " Radial"), CI(svmRResults$Sensitivity, ci=0.95)[2], CI(svmRResults$Specificity, ci=0.95)[2], CI(svmRResults$AccuracyPValue, ci=0.95)[2],CI(svmRResults$Accuracy, ci=0.95)[2])))
  output <- rbind(output, as.data.frame(  cbind(paste(GrfgDefinition, " rf"), CI(rfResults$Sensitivity, ci=0.95)[2], CI(rfResults$Specificity, ci=0.95)[2], CI(rfResults$AccuracyPValue, ci=0.95)[2],CI(rfResults$Accuracy, ci=0.95)[2])))
  
  colnames(output) <- c("Definition","Sensitivity", "Specificity", "PValue","Accuracy" )
  
 
  par(mfrow = c(1, 1))

  boxplot(nnetResults$Sensitivity, nnetResults$Specificity ,
          svmLResults$Sensitivity, svmLResults$Specificity,
          svmPResults$Sensitivity, svmPResults$Specificity ,
          svmRResults$Sensitivity, svmRResults$Specificity,
          rfResults$Sensitivity, rfResults$Specificity,
          main=paste(" Model Comparision for" , GrfgDefinition),
          at= c(1,2,4,5,7,8,10,11,13,14),
          names=c(" ",  "nnet",
                  " ", "svmLinear",
                  " ",  "svmPoly",
                  " ",  "svmRadial",
                  " ",  "rf"
                 ),
          col = c("orange","red"),
          las=2,
          border = "brown",
          horizontal = TRUE,
          notch = FALSE)
  mtext(paste("Orange:Sensitivity   Red: Specificity"))

   
  output 
}

CalculateAccuracyResultsv2 <- function(accuracy_result, trail_code,GrfgDefinition)
{
  # models <- c( "svmLinear","ranger", "svmPoly","nnet","naive_bayes","svmRadialGrid","adaBoost")
  output_sub<-data.frame(V1<-character(0),V2<-numeric(0) ,  V3<-numeric(0) , V4<-numeric(0) ,      V5<-numeric(0) , V6<-numeric(0) ,  V7 <-numeric(0) )
  
  
  accuracy_result <- accuracy_result$result[accuracy_result$result$trialCode %in% trail_code,]
  
  cat(paste("nnet \n"))
  nnetResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "nnet", ]  
  
  cat(paste("SVM Linear \n"))
  svmLResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "svmLinear", ]  
  
  cat(paste("SVM Poly \n"))
  svmPResults <-accuracy_result[accuracy_result$`model_list[j]` %in% "svmPoly", ]  
  
  cat(paste("SVM radial-1 \n"))
  svmRResults <- accuracy_result[accuracy_result$`model_list[j]` %in% c("svmRadialGrid","svmRadial"), ]  
  
  cat(paste("svmPResultsSVM radial-2 \n"))
  svmRResults2 <- accuracy_result[accuracy_result$`model_list[j]` %in% c("svmRadialGrid3"), ] 
  
  cat(paste("Ranger \n")) 
  rfResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "ranger", ] 
  
  cat(paste("Naive Bayes \n")) 
  nbResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "naive_bayes", ]  
  
  cat(paste("Adaboost \n")) 
  adaResults <- accuracy_result[accuracy_result$`model_list[j]` %in% "adaBoost", ] 
  
  cat(paste("Adaboost 2 \n")) 
  adaResults2 <- accuracy_result[accuracy_result$`model_list[j]` %in% "adaBoost2", ] 
  cat(paste("1")) 
  
  
  if(nrow(nnetResults[1])>0 ) { output_sub <-                    as.data.frame(  cbind(paste(GrfgDefinition, " nnet"),        CI(nnetResults$Sensitivity, ci=0.95)[2],  CI(nnetResults$Specificity, ci=0.95)[2],  CI(nnetResults$AccuracyPValue, ci=0.95)[2], CI(nnetResults$Accuracy, ci=0.95)[2],  CI(nnetResults$F1, ci=0.95)[2],   CI(nnetResults$Kappa, ci=0.95)[2] ))}
  cat(paste("2")) 
  if(nrow(svmLResults[1])>0 ) { output_sub <-  rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " Linear"),      CI(svmLResults$Sensitivity, ci=0.95)[2],  CI(svmLResults$Specificity, ci=0.95)[2],  CI(svmLResults$AccuracyPValue, ci=0.95)[2], CI(svmLResults$Accuracy, ci=0.95)[2],  CI(svmLResults$F1, ci=0.95)[2],  CI(svmLResults$Kappa, ci=0.95)[2] )))}
  cat(paste("3")) 
  if(nrow(svmPResults[1])>0 ) { output_sub <-  rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " Poly"),        CI(svmPResults$Sensitivity, ci=0.95)[2],  CI(svmPResults$Specificity, ci=0.95)[2],  CI(svmPResults$AccuracyPValue, ci=0.95)[2], CI(svmPResults$Accuracy, ci=0.95)[2],  CI(svmPResults$F1, ci=0.95)[2],  CI(svmPResults$Kappa, ci=0.95)[2] )))}
  cat(paste("4")) 
  if(nrow(svmRResults[1])>0 ) { output_sub <- rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " Radial"),      CI(svmRResults$Sensitivity, ci=0.95)[2],  CI(svmRResults$Specificity, ci=0.95)[2],  CI(svmRResults$AccuracyPValue, ci=0.95)[2], CI(svmRResults$Accuracy, ci=0.95)[2],  CI(svmRResults$F1, ci=0.95)[2],  CI(svmRResults$Kappa, ci=0.95)[2] )))}
  cat(paste("5")) 
  if(nrow(svmRResults2[1])>0 ) { output_sub <- rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " Radial3"),     CI(svmRResults2$Sensitivity, ci=0.95)[2], CI(svmRResults2$Specificity, ci=0.95)[2],CI(svmRResults2$AccuracyPValue, ci=0.95)[2], CI(svmRResults2$Accuracy, ci=0.95)[2], CI(svmRResults2$F1, ci=0.95)[2], CI(svmRResults2$Kappa, ci=0.95)[2] )))}
  cat(paste("6")) 
  if(nrow(rfResults[1])>0 ) { output_sub <-    rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " ranger"),      CI(rfResults$Sensitivity, ci=0.95)[2],    CI(rfResults$Specificity, ci=0.95)[2],    CI(rfResults$AccuracyPValue, ci=0.95)[2],   CI(rfResults$Accuracy, ci=0.95)[2],    CI(rfResults$F1, ci=0.95)[2],    CI(rfResults$Kappa, ci=0.95)[2] )))}
  cat(paste("7")) 
  if(nrow(nbResults[1])>0 ) {  output_sub <-   rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " naive_bayes"), CI(nbResults$Sensitivity, ci=0.95)[2],    CI(nbResults$Specificity, ci=0.95)[2],    CI(nbResults$AccuracyPValue, ci=0.95)[2],   CI(nbResults$Accuracy, ci=0.95)[2],    CI(nbResults$F1, ci=0.95)[2],    CI(nbResults$Kappa, ci=0.95)[2] )))}
  cat(paste("8")) 
   if(nrow(adaResults[1])>0 ) {  output_sub <- rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " adaBoost"),    CI(adaResults$Sensitivity, ci=0.95)[2],   CI(adaResults$Specificity, ci=0.95)[2],   CI(adaResults$AccuracyPValue, ci=0.95)[2],  CI(adaResults$Accuracy, ci=0.95)[2],   CI(adaResults$F1, ci=0.95)[2],   CI(adaResults$Kappa, ci=0.95)[2] )))}
  cat(paste("9")) 
  if(nrow(adaResults2[1])>0 ) { output_sub <-  rbind(output_sub, as.data.frame(  cbind(paste(GrfgDefinition, " adaBoost2"),   CI(adaResults2$Sensitivity, ci=0.95)[2],  CI(adaResults2$Specificity, ci=0.95)[2],  CI(adaResults2$AccuracyPValue, ci=0.95)[2], CI(adaResults2$Accuracy, ci=0.95)[2],  CI(adaResults2$F1, ci=0.95)[2],  CI(adaResults2$Kappa, ci=0.95)[2] )))} 
  
  colnames(output_sub) <- c("Definition","Sensitivity", "Specificity", "PValue","Accuracy",'F Value', 'Kappa' )
  
  
  par(mfrow = c(1, 1))

  boxplot(nnetResults$Sensitivity, nnetResults$Specificity ,
          svmLResults$Sensitivity, svmLResults$Specificity,
          svmPResults$Sensitivity, svmPResults$Specificity ,
          svmRResults$Sensitivity, svmRResults$Specificity,
          svmRResults2$Sensitivity, svmRResults2$Specificity,
          rfResults$Sensitivity, rfResults$Specificity,
          nbResults$Sensitivity, nbResults$Specificity,
          adaResults$Sensitivity, adaResults$Specificity,
          adaResults2$Sensitivity, adaResults2$Specificity,
          main=paste(" Model Comparision for" , GrfgDefinition),
          at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),
          names=c(" ",  "nnet",
                  " ", "svmLinear",
                  " ",  "svmPoly",
                  " ",  "svmRadial",
                  " ",  "svmRadial2",
                  " ",  "ranger",
                  " ",  "naive_bayes",
                  " ",  "adaBoost",
                  " ",  "adaBoost2"
          ),
          col = c("orange","red"),
          las=2,
          border = "brown",
          horizontal = TRUE,
          notch = FALSE)
  mtext(paste("Orange:Sensitivity   Red: Specificity"))

  
  output_sub 
}


CompareAccuracyResults <- function(results,CompareDataFrame,Definition)
{
  outpt <- NULL
  ci <- 0.95
  # CompareDataFrame$trail_code, CompareDataFrame$model_type, CompareDataFrame$
 
  modelResultsList <- NULL

  for(i in 1:nrow(CompareDataFrame)){
    #cat(paste(i,"-item \n"))
    accuracy_result <- results[results$trialCode %in% CompareDataFrame$trail_code[i] ,]
    head(accuracy_result)
    modelResults <- accuracy_result[accuracy_result$`model_list[j]` %in% paste(CompareDataFrame$model_type[i],"",sep="") , ]  
   
    if (nrow(modelResults)==0)
    {
      modelResults <- accuracy_result[accuracy_result$`model_list[j]` %in% CompareDataFrame$model_type[i] , ]  
      
    }
    
  
    head(accuracy_result)
    df <- as.data.frame(cbind(
      paste(Definition, CompareDataFrame$model_type[i] ),
      CI(modelResults$Sensitivity, ci = ci)[2],
      CI(modelResults$Specificity, ci = ci)[2],
      CI(modelResults$Accuracy, ci = ci)[2],
      CI(modelResults$AccuracyPValue, ci = ci)[2],
      CI(modelResults$Kappa, ci = ci)[2],
      CI(modelResults$F1, ci = ci)[2],
      CI(modelResults$McnemarPValue, ci = ci)[2]
    ))
    
    
    head(modelResults)
    if(i ==1){
      outpt <- df
      modelResultsList <-  cbind(i,modelResults)
    } else { 
      outpt <- rbind(outpt, df) 
      modelResultsList <- rbind(modelResultsList,cbind(i,modelResults))
      }
    
   # cat(paste("Completed \n"))
   
  }
 

  head(modelResultsList)
  
  colnames(outpt) <- c("Definition","Sensitivity", "Specificity", "Accuracy" , "PValue", "Kappa", "FScore","McnemarPValue")
  
  
  boxplot(modelResultsList[modelResultsList$i %in% 1,]$Sensitivity, modelResultsList[modelResultsList$i %in% 1,]$Specificity ,modelResultsList[modelResultsList$i %in% 1,]$Accuracy,modelResultsList[modelResultsList$i %in% 1,]$F1,
          modelResultsList[modelResultsList$i %in% 2,]$Sensitivity, modelResultsList[modelResultsList$i %in% 2,]$Specificity,modelResultsList[modelResultsList$i %in% 2,]$Accuracy,modelResultsList[modelResultsList$i %in% 2,]$F1,
           modelResultsList[modelResultsList$i %in% 3,]$Sensitivity, modelResultsList[modelResultsList$i %in% 3,]$Specificity ,modelResultsList[modelResultsList$i %in% 3,]$Accuracy,modelResultsList[modelResultsList$i %in% 3,]$F1,
           modelResultsList[modelResultsList$i %in% 4,]$Sensitivity, modelResultsList[modelResultsList$i %in% 4,]$Specificity,modelResultsList[modelResultsList$i %in% 4,]$Accuracy,modelResultsList[modelResultsList$i %in% 4,]$F1,
          modelResultsList[modelResultsList$i %in% 5,]$Sensitivity, modelResultsList[modelResultsList$i %in% 5,]$Specificity ,modelResultsList[modelResultsList$i %in% 5,]$Accuracy,modelResultsList[modelResultsList$i %in% 5,]$F1,
           modelResultsList[modelResultsList$i %in% 6,]$Sensitivity, modelResultsList[modelResultsList$i %in% 6,]$Specificity ,modelResultsList[modelResultsList$i %in% 6,]$Accuracy,modelResultsList[modelResultsList$i %in% 6,]$F1,
          main= Definition,
          at= c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29),
          names=c(" "," ", CompareDataFrame$def[1]," ",
                  " "," ", CompareDataFrame$def[2]," "
                  , " "," ", CompareDataFrame$def[3]," ",
                   " "," ", CompareDataFrame$def[4]," ",
                   " "," ", CompareDataFrame$def[5]," ",
                   " "," ", CompareDataFrame$def[6]," "
          ),
          col = rainbow(4, alpha=0.2), #c("grey", rgb(0.1,0.1,0.7,0.5)),
          las=2,
          border = rainbow(4, v=0.4),
          horizontal = FALSE,
          notch = FALSE,
          ylim = c(0, 1),
          xlab="Red:Sensitivity  Green:Specificity  Blue:Accuracy  Purple:FScore " )
 
  
  
  outpt 
}
CompareAccuracyResult_9Set <- function(results,CompareDataFrame,Definition)
{
  outpt <- NULL
  ci <- 0.95
  # CompareDataFrame$trail_code, CompareDataFrame$model_type, CompareDataFrame$
  
  modelResultsList <- NULL
  
  for(i in 1:nrow(CompareDataFrame)){
    #cat(paste(i,"-item \n"))
    accuracy_result <- results[results$trialCode %in% CompareDataFrame$trail_code[i] ,]
    head(accuracy_result)
    modelResults <- accuracy_result[accuracy_result$`model_list[j]` %in% paste(CompareDataFrame$model_type[i],"",sep="") , ]  
    
    if (nrow(modelResults)==0)
    {
      modelResults <- accuracy_result[accuracy_result$`model_list[j]` %in% CompareDataFrame$model_type[i] , ]  
      
    }
    
    
    head(accuracy_result)
    df <- as.data.frame(cbind(
      paste(Definition, CompareDataFrame$model_type[i] ),
      CI(modelResults$Sensitivity, ci = ci)[2],
      CI(modelResults$Specificity, ci = ci)[2],
      CI(modelResults$Accuracy, ci = ci)[2],
      CI(modelResults$AccuracyPValue, ci = ci)[2],
      CI(modelResults$Kappa, ci = ci)[2],
      CI(modelResults$F1, ci = ci)[2],
      CI(modelResults$McnemarPValue, ci = ci)[2]
    ))
    
    
    head(modelResults)
    if(i ==1){
      outpt <- df
      modelResultsList <-  cbind(i,modelResults)
    } else { 
      outpt <- rbind(outpt, df) 
      modelResultsList <- rbind(modelResultsList,cbind(i,modelResults))
    }
    
    # cat(paste("Completed \n"))
    
  }
  
  
  head(modelResultsList)
  
  colnames(outpt) <- c("Definition","Sensitivity", "Specificity", "Accuracy" , "PValue", "Kappa", "FScore","McnemarPValue")
  
  
  boxplot(modelResultsList[modelResultsList$i %in% 1,]$Sensitivity, modelResultsList[modelResultsList$i %in% 1,]$Specificity ,modelResultsList[modelResultsList$i %in% 1,]$Accuracy,modelResultsList[modelResultsList$i %in% 1,]$F1,
          modelResultsList[modelResultsList$i %in% 2,]$Sensitivity, modelResultsList[modelResultsList$i %in% 2,]$Specificity,modelResultsList[modelResultsList$i %in% 2,]$Accuracy,modelResultsList[modelResultsList$i %in% 2,]$F1,
          modelResultsList[modelResultsList$i %in% 3,]$Sensitivity, modelResultsList[modelResultsList$i %in% 3,]$Specificity ,modelResultsList[modelResultsList$i %in% 3,]$Accuracy,modelResultsList[modelResultsList$i %in% 3,]$F1,
          modelResultsList[modelResultsList$i %in% 4,]$Sensitivity, modelResultsList[modelResultsList$i %in% 4,]$Specificity,modelResultsList[modelResultsList$i %in% 4,]$Accuracy,modelResultsList[modelResultsList$i %in% 4,]$F1,
          modelResultsList[modelResultsList$i %in% 5,]$Sensitivity, modelResultsList[modelResultsList$i %in% 5,]$Specificity ,modelResultsList[modelResultsList$i %in% 5,]$Accuracy,modelResultsList[modelResultsList$i %in% 5,]$F1,
          modelResultsList[modelResultsList$i %in% 6,]$Sensitivity, modelResultsList[modelResultsList$i %in% 6,]$Specificity ,modelResultsList[modelResultsList$i %in% 6,]$Accuracy,modelResultsList[modelResultsList$i %in% 6,]$F1,
          modelResultsList[modelResultsList$i %in% 7,]$Sensitivity, modelResultsList[modelResultsList$i %in% 7,]$Specificity,modelResultsList[modelResultsList$i %in% 7,]$Accuracy,modelResultsList[modelResultsList$i %in% 7,]$F1,
          modelResultsList[modelResultsList$i %in% 8,]$Sensitivity, modelResultsList[modelResultsList$i %in% 8,]$Specificity ,modelResultsList[modelResultsList$i %in% 8,]$Accuracy,modelResultsList[modelResultsList$i %in% 8,]$F1,
          modelResultsList[modelResultsList$i %in% 9,]$Sensitivity, modelResultsList[modelResultsList$i %in% 9,]$Specificity ,modelResultsList[modelResultsList$i %in% 9,]$Accuracy,modelResultsList[modelResultsList$i %in% 9,]$F1,
          main= Definition,
          at= c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29,31,32,33,34,36,37,38,39,41,42,43,44),
          names=c(" "," ", CompareDataFrame$def[1]," ",
                  " "," ", CompareDataFrame$def[2]," " , 
                  " "," ", CompareDataFrame$def[3]," ",
                  " "," ", CompareDataFrame$def[4]," ",
                  " "," ", CompareDataFrame$def[5]," ",
                  " "," ", CompareDataFrame$def[6]," ",
                  " "," ", CompareDataFrame$def[7]," ",
                  " "," ", CompareDataFrame$def[8]," ",
                  " "," ", CompareDataFrame$def[9]," "
          ),
          col = rainbow(4, alpha=0.2), #c("grey", rgb(0.1,0.1,0.7,0.5)),
          las=2,
          border = rainbow(4, v=0.4),
          horizontal = FALSE,
          notch = FALSE,
          ylim = c(0, 1),
          xlab="1:miRNA  2:mRNA  3:Methylation  4:miRNA & mRNA  5:mRNA & Methylation   6:miRNA & Methylation   7:miRNA mRNA & Methylation   8:miRNA mRNA & Low Methylation    9:miRNA mRNA & H Methylation" )
  
  
  
  outpt 
}
tunnePlot <-function(tunneResults) {
  
  boxplot(tunneResults[tunneResults$size %in% 1,]$Accuracy,
          tunneResults[tunneResults$size %in% 2,]$Accuracy,
          tunneResults[tunneResults$size %in% 3,]$Accuracy,
          tunneResults[tunneResults$size %in% 4,]$Accuracy,
          tunneResults[tunneResults$size %in% 5,]$Accuracy,
          tunneResults[tunneResults$size %in% 6,]$Accuracy,
          tunneResults[tunneResults$size %in% 7,]$Accuracy,
          tunneResults[tunneResults$size %in% 8,]$Accuracy,
          tunneResults[tunneResults$size %in% 9,]$Accuracy,
          tunneResults[tunneResults$size %in% 10,]$Accuracy,
          main= "Tunning for NNET-Size",
          at= c(1,2,3,4,5,6,7,8,9,10),
          names=c(1:10),
          col = rainbow(1, alpha=0.2), #c("grey", rgb(0.1,0.1,0.7,0.5)),
          las=2,
          border = rainbow(1, v=0.4),
          horizontal = FALSE,
          notch = FALSE,
          ylim = c(0, 1) )
  
  # 0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.04,0.06,0.08,0.1,1,3
  boxplot(tunneResults[tunneResults$decay %in% 0.001,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.002,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.004,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.006,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.008,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.01,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.02,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.04,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.06,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.08,]$Accuracy,
          tunneResults[tunneResults$decay %in% 0.1,]$Accuracy,
          tunneResults[tunneResults$decay %in% 1.00,]$Accuracy,
          tunneResults[tunneResults$decay %in% 3.00,]$Accuracy,
          tunneResults[tunneResults$decay %in% 5.00,]$Accuracy,
          main= "Tunning for NNET-Decay",
          at= c(1:14),
          names=c(" 0.001","0.002","0.004","0.006","0.008","0.01","0.02","0.04","0.06","0.08","0.1","1","3","5"),
          col = rainbow(1, alpha=0.2), #c("grey", rgb(0.1,0.1,0.7,0.5)),
          las=2,
          border = rainbow(1, v=0.4),
          horizontal = FALSE,
          notch = FALSE,
          ylim = c(0, 1) )
  
}


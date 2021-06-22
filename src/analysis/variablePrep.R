modelVariableSelecttions <- NULL
runn<- 1

# miRNA....

lvariables <- rbind(pool$miRNAFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables$probName)  
View(lvariablesList)

modelInpts <-  fxGetModelDataInputs (dataSample=dataSampleList[[runn]],mode="rf",variableList=lvariablesList,runtypes=c("sm","us"),usePrev=FALSE,data=data,preruuid=NULL)


modelVariableSelecttions$miRNA <- colnames(modelInpts$us$test)

#mRNA

lvariables <- rbind(pool$mRNAFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:500,]$probName)  
modelVariableSelecttions$mRNA <- lvariablesList 

#Methylation

lvariables <- rbind(pool$methylationFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
modelVariableSelecttions$Methy <- lvariablesList 

#mRNA and miRNA
lvariables <- rbind(pool$mRNAFeatures,pool$miRNAFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
modelVariableSelecttions$miRNA_mRNA <- lvariablesList 

#mRNA and Methylation
lvariables <- rbind(pool$mRNAFeatures, pool$methylationFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables$probName)  

modelInpts <-  fxGetModelDataInputs (dataSample=dataSampleList[[runn]],mode="rf",variableList=lvariablesList,runtypes=c("sm","us"),usePrev=FALSE,data=data,preruuid=NULL)
modelVariableSelecttions$mRNA_Methy<- colnames(modelInpts$sm$test)

#miRNA and Methylation
lvariables <- rbind(pool$miRNAFeatures, pool$methylationFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:500,]$probName)  
modelVariableSelecttions$miRNA_Methy <- lvariablesList 

# miRNA , mRNA and Methylation
lvariables <- rbind(pool$miRNAFeatures, pool$methylationFeatures)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
modelVariableSelecttions$miRNA_mRNA_Methy <- lvariablesList 


# miRNA , mRNA and Methylation L
lvariables <- rbind(pool$miRNAFeatures, pool$methylationFeatures.L)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
modelVariableSelecttions$miRNA_mRNA_Methy_L <- lvariablesList 

# miRNA , mRNA and Methylation L
lvariables <- rbind(pool$miRNAFeatures, pool$methylationFeatures.H)
setorder(lvariables,P_Val) 
lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
modelVariableSelecttions$miRNA_mRNA_Methy_H <- lvariablesList 


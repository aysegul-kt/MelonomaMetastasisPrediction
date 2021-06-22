library(readr)



Methylation.Path <- "~/Documents/Melonoma_TCGA-SKCM/DATASET/DNA_Methylation_Data"

datafiles<-list.files(path =Methylation.Path , full.names = TRUE, recursive = TRUE)

View(datafiles)



Methylation.metadata <- read_csv("~/Documents/Melonoma_TCGA-SKCM/DATASET/DNA_Methylation_Data/tcga_grch38_20200303_1.csv")

Methylation.data <- NULL

length(datafiles)
for (datafile in datafiles) {
    current_filename <- unlist(strsplit(datafile, "/"))[9]
    caseuuid <- Methylation.metadata[Methylation.metadata$file %in% current_filename, ]$case
    cat(nrow(Methylation.data)+1)
    tempdata <-
      read_delim(
        datafile,
        "\t",
        escape_double = FALSE,
        col_types = cols(
          CGI_Coordinate = col_skip(),
          Chromosome = col_skip(),
          End = col_skip(),
          Feature_Type = col_skip(),
          Gene_Type = col_skip(),
          Position_to_TSS = col_skip(),
          Start = col_skip(),
          Transcript_ID = col_skip()
        ),
        trim_ws = TRUE
      )
    
    item<-as.data.frame(  t(as.vector( tempdata$Beta_value)))
    colnames(item)<-as.character(tempdata$`Composite Element REF`)
    item$case_uuid <- c(caseuuid)
    
  if (is.null( Methylation.data ) ) { 
    Methylation.data <-  item 
  }else{
    Methylation.data <- rbind(Methylation.data,item)
  }
}




Methylation.data_cleaned  <- Filter(function(x) !(all(x=="")), Methylation.data)



colnames(Methylation.data_cleaned[ncol(Methylation.data_cleaned)-1])
colnames(Methylation.data_cleaned[1])

library(readr)

melonomaCases <- read_csv("data/raw/melonomaData.csv", 
                          col_types = cols_only(case_uuid = col_guess(), 
                                                sample_type = col_guess(),
                                                case_id = col_guess()))

MethylationData.Final <- merge(x = Methylation.data_cleaned, y = melonomaCases, by = "case_uuid", all = TRUE)

#logToFile("MethylationService","MethylationData.Final.csv",MethylationData.Final,"D")
save(MethylationData.Final,file="data/processed/MethylationData.Raw.Final.RData")
rm(Methylation.data_cleaned,Methylation.data_cleaned_norm, Methylation.data,Methylation.Path,datafiles,Methylation.metadata,caseuuid,current_filename,Methylation.data,Methylation.data_cleaned,melonomaCases)
##########
# region based ttest  has been removed
# Methylation.class1<-MethylationData.Final[MethylationData.Final$sample_type %in% "Metastatic", ]  
# Methylation.class2<-MethylationData.Final[MethylationData.Final$sample_type %in% "Primary Tumor", ]
# 
# Methylation.class1$sample_type<-rep(c(1),nrow(Methylation.class1))
# Methylation.class2$sample_type<-rep(c(0),nrow(Methylation.class2))
# 
# # 1: metastastic 0: Primary tumor
# Methylation.class1$sample_type<-factor(Methylation.class1$sample_type, levels = c(0, 1)) 
# Methylation.class2$sample_type<-factor(Methylation.class2$sample_type, levels = c(0, 1)) 
# Methylation.class1 <- as.data.frame(Methylation.class1)
# Methylation.class2 <- as.data.frame(Methylation.class2)
# rm(Methylation.class2,Methylation.class1)
# 
# Methylation.ttest.results<-significancyTest("ttest",Methylation.class1[,1:(ncol(Methylation.class1)-1)],Methylation.class2[1:(ncol(Methylation.class2)-1)])
# 
# View(Methylation.ttest.results)
# logToFile("MethylationService","Methylation.ttest.results.csv",Methylation.ttest.results,"L")
#

# for (i in 1:nrow(Methylation.ttest.results)) {
#   
#   if (Methylation.ttest.results$Diff[i]==1)  {
#     cat(paste("\n",i , colnames(Methylation.class1[i+1]), " ", Methylation.ttest.results$probName[i]))
#     if (colnames(Methylation.class1[i+1]) == Methylation.ttest.results$probName[i]) {
#       propResult<-TTest(colnames(Methylation.class1[i+1]), Methylation.class1[i+1],  Methylation.class2[i+1])
#       
#       Methylation.ttest.results$P_Val_Less[i]  <- propResult$P_Val_Less
#       Methylation.ttest.results$P_Val_Greater[i]  <- propResult$P_Val_Greater
#       cat(paste(" L: ", propResult$P_Val_Less, " G: ",propResult$P_Val_Greater ))
#    }
#     }
#   
# }

# logToFile("MethylationService","Methylation.ttest.results.csv",Methylation.ttest.results,"L")
#Degree0001 <- Methylation.ttest.results #subset(Methylation.ttest.results, P_Val<0.001)


##############

### Sub set that will be used for the rest of the study
MethylationSelection <- MethylationData.Final #subset(MethylationData.Final,select = as.vector( Degree0001$probName ))  # all  data will be used
MethylationSelection$class <- MethylationData.Final$sample_type
MethylationSelection$case_uuid <-  MethylationData.Final$case_uuid
MethylationSelection$case_id <-  MethylationData.Final$case_id
### We will use one of the data sample as template to gerenare Element Region Ref to Gene Symbol...

MappingTemplate <-
  read_delim(
    datafiles[1],
    "\t",
    escape_double = FALSE,
    col_types = cols(
      CGI_Coordinate = col_skip(),
      Chromosome = col_skip(),
      End = col_skip(),
      Feature_Type = col_skip(),
      Gene_Type = col_skip(),
      Position_to_TSS = col_skip(),
      Start = col_skip(),
      Transcript_ID = col_skip(),
      Beta_value = col_skip()
    ),
    trim_ws = TRUE
  )

MappingTemplate <- MappingTemplate[MappingTemplate$`Composite Element REF` %in% Degree0001$probName,]
save(MappingTemplate,file="data/processed/Methylation.MappingTemplate.RData")

summary(MappingTemplate)
View(MappingTemplate)

MethylationMapping <-  data.frame(CEF =rep(NA,nrow(MappingTemplate)*4) ,hgnc_symbol =rep(NA,nrow(MappingTemplate)*4)  )  
y <- 1                   

for(i in 1:nrow(MappingTemplate) ){ # 
  cat(i)  
  geneSymbols<- unlist(strsplit(MappingTemplate$Gene_Symbol[i], ";"))
  geneSymbols <- geneSymbols[!duplicated(geneSymbols)]
  for (variable in geneSymbols) {
    if (variable != ".")
    {
      MethylationMapping$CEF[y] <-MappingTemplate$`Composite Element REF`[i]
      MethylationMapping$hgnc_symbol[y] <-variable  
      y <- y +1
    }
   
  }
  cat("\014")

}

summary(MethylationMapping)
View(MethylationMapping)

save(MethylationMapping,file="data/processed/Methylation.Mapping.RData")
# summation of beta values for same gene on different regions....

mylist <- unique(MethylationMapping$hgnc_symbol)

geneMethylationDataSet<-  data.frame(case_uuid =rep(NA,nrow(MethylationSelection)) ,class =rep(NA,nrow(MethylationSelection))  )  
geneMethylationDataSet$case_uuid <- MethylationSelection$case_uuid 
geneMethylationDataSet$class <- MethylationSelection$class

columnNameList <- colnames(geneMethylationDataSet)
y<-3
# size 15450

for (variable in mylist) {
 
  cat(y-2)
  temp <- MethylationMapping[MethylationMapping$hgnc_symbol %in% variable,]
  tempSubSet <- subset(MethylationData.Final,select = as.vector( temp$CEF ))
  geneMethylationDataSet$newColumn <- rowMeans(tempSubSet,na.rm = FALSE, dims = 1) 
  columnNameList[y] <- variable
  colnames(geneMethylationDataSet) <-  columnNameList
  y <- y+1
  cat("\014")
  
}
y
geneMethylationDataSet <- geneMethylationDataSet[,1:34018]  
geneMethylationDataSet$case_id <- MethylationSelection$case_id 
geneMethylationDataSet_norm <- dataNormalization(3,geneMethylationDataSet, ncol(geneMethylationDataSet)-1 ) # last item is case uuid

save(geneMethylationDataSet_norm,file="data/processed/MethylationData.Gene.RData")

rm(MethylationSelection,geneMethylationDataSet,MethylationMapping,MappingTemplate)
#logToFile("MethylationService","geneMethylationDataSet.csv",geneMethylationDataSet_norm,"L")

GeneMethylation.class1<-geneMethylationDataSet_norm[geneMethylationDataSet_norm$class %in% "Metastatic", ]  
GeneMethylation.class2<-geneMethylationDataSet_norm[geneMethylationDataSet_norm$class %in% "Primary Tumor", ]

GeneMethylation.class1$sample_type<-rep(c(1),nrow(GeneMethylation.class1))
GeneMethylation.class2$sample_type<-rep(c(0),nrow(GeneMethylation.class2))

# 1: metastastic 0: Primary tumor
GeneMethylation.class1$sample_type<-factor(GeneMethylation.class1$sample_type, levels = c(0, 1)) 
GeneMethylation.class2$sample_type<-factor(GeneMethylation.class2$sample_type, levels = c(0, 1)) 
GeneMethylation.class1 <- as.data.frame(GeneMethylation.class1)
GeneMethylation.class2 <- as.data.frame(GeneMethylation.class2)


GeneMethylation.ttest.results<-significancyTest("ttest",GeneMethylation.class1[,3:(ncol(GeneMethylation.class1)-1)],GeneMethylation.class2[,3:(ncol(GeneMethylation.class2)-1)])
GeneMethylation.ttest.results <- na.omit( GeneMethylation.ttest.results)
View(GeneMethylation.ttest.results)

logToFile("MethylationService","GeneMethylation.ttest.results.csv",GeneMethylation.ttest.results,"L")
save(GeneMethylation.ttest.results,file="data/processed/GeneMethylation.ttest.results.RData")

rm(GeneMethylation.class1,GeneMethylation.class2)
rm(Degree0001)

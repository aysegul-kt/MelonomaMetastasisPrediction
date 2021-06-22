# miRNA- mRNA dataset Prepare

library("readr")

# setwd("~/Documents/Melonoma_TCGA-SKCM")

manifest <- read_csv("~/Documents/Melonoma_TCGA-SKCM/DATASET/manifest.csv")

cases <- read_csv("~/Documents/Melonoma_TCGA-SKCM//DATASET/cases.csv")


for (j in 1:470) 
{
  
  caseItem <- as.character(cases$CaseId[j])
  dirs<-list.dirs(path = paste("./DATASET/data/", caseItem,sep=""), full.names = TRUE, recursive = FALSE)
  
  k<-1
  for (k in 1:length(dirs))
  {
    i<-1
    item<- gsub(paste("./DATASET/data/", caseItem,"/",sep=""),"",dirs[k])
    
    for (i in 1:3260)
    {
      id <-  as.character( manifest$id[i])
      
      if(id==item)
      {
        manifest$caseId[i]<- caseItem
        manifest$path[i]<-paste("./DATASET/data/", caseItem,"/",sep="")
      }
    }
  }
}

write.csv(file="./2018_DEC/melonomaData_fileset.csv",manifest)

melonomaData_fileset <- read_csv("~/Documents/Melonoma_TCGA-SKCM/2018_DEC/melonomaData_fileset.csv")
mirna<-0
rnafpkm<-0
rnafpkmuq<-0
hseq<-0

for (j in 1:3260) 
{

filename<- as.character( melonomaData_fileset$filename[j])
case_id<- as.character( melonomaData_fileset$caseId[j])
path<- as.character( melonomaData_fileset$path[j])
id<- as.character( melonomaData_fileset$id[j])

metaDataRow<-cbind(case_id,path,id,filename)

#mirna 
if(length(grep( "mirnas.quantification", filename)) >0){
  mirna<-mirna+1
  dataFile<-read.delim(paste(path,id,"/",filename,sep=""))
  miRNA_Data<-t(as.vector( dataFile$reads_per_million_miRNA_mapped))
  colnames(miRNA_Data)<-as.character( dataFile$miRNA_ID)
  
  #bind
  d<-cbind(metaDataRow,miRNA_Data)
  #as.vector( myMetaData[2,])
  #as.vector(miRNA_Data)
  if(mirna>1)
  {
    miRNA_Data_Segment<-rbind(miRNA_Data_Segment,d)
  }else{
    miRNA_Data_Segment<-d
  }
}


#htseq.counts

if(length(grep( "htseq.counts", filename)) >0){
 
  hseq<-hseq+1
  dataFile<-read.delim(paste(path,id,"/",filename,sep=""),header=FALSE)
  dataRowNames<- as.character(dataFile$V1)
  htseqDATA<-t(as.vector( dataFile[2]))
  colnames(htseqDATA)<-dataRowNames
  
  #bind
  d<-cbind(metaDataRow,htseqDATA)
  #as.vector( myMetaData[2,])
  #as.vector(miRNA_Data)
  if(hseq>1)
  {
    htseqDATA_Segment<-rbind(htseqDATA_Segment,d)
  }else{
    htseqDATA_Segment<-d
  }
} 

}
# 
# #FPKM-UQ.txt
# if(length(grep( "FPKM-UQ.txt", filename)) >0){
#   rnafpkmuq<-rnafpkmuq+1
#   dataFile<-read.delim(paste(path,id,"/",filename,sep=""),header=FALSE)
#   dataRowNames<- as.character(dataFile$V1)
#   RNA_FPKM_UQ_Data<-t(as.vector( dataFile[2]))
#   colnames(RNA_FPKM_UQ_Data)<-dataRowNames
#   
#   #bind
#   d<-cbind(metaDataRow,RNA_FPKM_UQ_Data)
#   #as.vector( myMetaData[2,])
#   #as.vector(miRNA_Data)
#   if(rnafpkmuq>2)
#   {
#     RNA_FPKM_UQ_Data_Segment<-rbind(RNA_FPKM_UQ_Data_Segment,d)
#   }else{
#     RNA_FPKM_UQ_Data_Segment<-d
#   }
# }
# }
# 
# #FPKM.txt
# 
# if(length(grep( "FPKM.txt", filename)) >0){
#   rnafpkm<-rnafpkm+1
#   dataFile<-read.delim(paste(path,id,"/",filename,sep=""),header=FALSE)
#   dataRowNames<- as.character(dataFile$V1)
#   RNA_FPKM_Data<-t(as.vector( dataFile[2]))
#   colnames(RNA_FPKM_Data)<-dataRowNames
#   
#   #bind
#   d<-cbind(metaDataRow,RNA_FPKM_Data)
#   #as.vector( myMetaData[2,])
#   #as.vector(miRNA_Data)
#   if(rnafpkm>2)
#   {
#     RNA_FPKM_Data_Segment<-rbind(RNA_FPKM_Data_Segment,d)
#   }else{
#     RNA_FPKM_Data_Segment<-d
#   }
#   


#isoforms.quantification.txt
#do not use..

# write.csv(file="2018_DEC/RNA_FPKM_UQ_Data_Segment.csv",RNA_FPKM_UQ_Data_Segment)
# write.csv(file="2018_DEC/htseqDATA_Segment.csv",htseqDATA_Segment)
# write.csv(file="2018_DEC/RNA_FPKM_Data_Segment.csv",RNA_FPKM_Data_Segment)
# write.csv(file="2018_DEC/miRNA_Data_Segment.csv",miRNA_Data_Segment)


miRNA_list <-  as.vector( colnames(miRNA_Data))
mRNA_list <-  as.vector( colnames(htseqDATA))


write.csv(file="./2018_DEC/miRNA_list.csv",miRNA_list)
write.csv(file="./2018_DEC/mRNA_list.csv",mRNA_list)

# melonomaData<-melonomaData[melonomaData$sample_type %in% c("Metastatic","Primary Tumor"),]
cat(paste("\n Clone miRNA List"))
miRNA_list_clone<-miRNA_list

for (i in 1:length(miRNA_list))
{
 if ( (nchar(miRNA_list[i])- nchar(gsub("-","",miRNA_list[i])))>2)
 {
   miRNA_list_clone[i]<- substr(miRNA_list[i], 1, nchar(miRNA_list[i])-2)
   cat(paste("..", i, ".."))
 }
  
}
#write.csv(file="2018_DEC/miRNA_list_clone.csv",miRNA_list_clone)


# find classes of datas
melonomaData <- read_csv("~/Documents/Melonoma_TCGA-SKCM/DATASET/melonomaData.csv")

miRNA_Data_Segment<-as.data.frame(miRNA_Data_Segment)
htseqDATA_Segment<-as.data.frame(htseqDATA_Segment)
clinicRow<- matrix(nrow=1, ncol=19)
colnames(clinicRow) <- colnames(melonomaData[,4:22])

# mirna
cat(paste("\n  Prepare: miRNA_Data_Segment Total:",nrow(miRNA_Data_Segment)))
for (i in 1:nrow(miRNA_Data_Segment))
{
  cat(paste(".", i, "."))

 r<- melonomaData[melonomaData$case_id%in%  miRNA_Data_Segment$case_id[i],4:22]
 samplecount<-nrow(miRNA_Data_Segment[miRNA_Data_Segment$case_id%in%  miRNA_Data_Segment$case_id[i],])
 if (samplecount==1 && nrow(r)>0)
  {
    myRow <- cbind(r,miRNA_Data_Segment[i,])
  } else {
    myRow <- cbind(clinicRow,miRNA_Data_Segment[i,])
  }
 
  if (i>1){
     miRNA_data_final <- rbind(miRNA_data_final,myRow)
  }else{
     miRNA_data_final <- myRow
   }
}

#rna
cat(paste("\n  Prepare: mRNA_Data_Segment: Total",nrow(htseqDATA_Segment)))
for (i in 1:nrow(htseqDATA_Segment))
{
  cat(paste("\n", i))
  r<- melonomaData[melonomaData$case_id%in%  htseqDATA_Segment$case_id[i],4:22]
  samplecount<-nrow(htseqDATA_Segment[htseqDATA_Segment$case_id%in%  htseqDATA_Segment$case_id[i],])
  
  if (samplecount==1 && nrow(r)>0)
  {
    myRow <- cbind(r,htseqDATA_Segment[i,])
  } else {
    myRow <- cbind(clinicRow,htseqDATA_Segment[i,])
  }
  
  if (i>1){
    colnames(myRow)<-colnames(htseq_data_final)
    htseq_data_final <- rbind(htseq_data_final,myRow)
  }else{
    htseq_data_final <- myRow
  }

}

write.csv(file="./2018_DEC/htseq_data_final.csv",htseq_data_final)
write.csv(file="./2018_DEC/miRNA_data_final.csv",miRNA_data_final)

# # remove null sample types....
# 
# # miRNA_data_step_one<-miRNA_data_final[miRNA_data_final$sample_type %in% c("Metastatic","Primary Tumor"), ]  
# # htseq_data_step_one<-htseq_data_final[htseq_data_final$sample_type %in% c("Metastatic","Primary Tumor"), ]  
# 
# miRNA_data_step_one<-removeNullClass(miRNA_data_final)
# htseq_data_step_one<-removeNullClass(htseq_data_final)
# 
# write.csv(file="./2018_DEC/miRNA_data_step_one.csv",miRNA_data_step_one)
# write.csv(file="./2018_DEC/dhtseq_data_step_one.csv",htseq_data_step_one)
# 
# 
# miRNA_data_step_two<-variableSelection(miRNA_data_step_one,"mirna")
# write.csv(file="./2018_DEC/miRNA_data_step_two.csv",miRNA_data_step_two)
# 
# # htseq_data_step_one$sample_type<-htseq_data_step_one[6]
# htseq_data_step_one$sample_type<-as.factor((unlist(htseq_data_step_one$sample_type)))
# htseq_data_step_two<-variableSelection(htseq_data_step_one,"mrna")
# write.csv(file="./2018_DEC/dhtseq_data_step_two.csv",htseq_data_step_two)


###############################################
# ADI              :  Data base operations 
# YAZAR            : AYSEGUL KUTLAY
# ACIKLAMA         :    
# GUNCELLEME TAR?H?: 06.06.2019
# VERSIYON         : 02
################################################

# Database 

miRTarBase_MTI <- read_csv("./mirTARbase/miRTarBase_MTI.csv")
View(miRTarBase_MTI)


miRTarBase_MTI$miRNA<-tolower(miRTarBase_MTI$miRNA)
# miRNA_list_clone<-as.data.frame(miRNA_list_clone)
# targetList<-miRTarBase_MTI[miRTarBase_MTI%in% miRNA_list_clone$miRNA_list_clone ,]
targetList<-miRTarBase_MTI[miRTarBase_MTI$`Species (miRNA)` %in% "Homo sapiens" ,]

logToFile("DatabaseOps","targetList.csv",targetList,"D")

# write.csv(file="targetList.csv",targetList)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("biomaRt"))

 # source("http://bioconductor.org/biocLite.R")
 # biocLite("biomaRt")
 library(biomaRt)

 
 
# Ensembl US West: http://uswest.ensembl.org/index.html
# Ensembl US East: http://useast.ensembl.org/index.html
# Ensembl Asia: http://asia.ensembl.org/index.html

listMarts(host="ensembl.org")
listEnsembl(host="useast.ensembl.org")
ensembl = useEnsembl(biomart="ensembl",host="ensembl.org")
listDatasets(ensembl)
human = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="ensembl.org") 
targetList$`Target Gene`
geneList<-unique(  as.vector( targetList$`Target Gene`))


for (i in 1:length(geneList))
{
  cat(i)
  hgnc_item <- getBM(attributes=c('ensembl_gene_id','ensembl_gene_id_version','ensembl_transcript_id','ensembl_transcript_id_version','hgnc_symbol','hgnc_id','description'),filters = 'hgnc_symbol', values = geneList, mart = human)
  if(i>2)
  {
    hgnc_list<-rbind(hgnc_list,hgnc_item)
  }else{
    hgnc_list<-hgnc_item
  }
  cat("\014") 
}


logToFile("DatabaseOps","hgnc_list.csv",hgnc_list,"D")

#write.csv(file="hgnc_list.csv",hgnc_list)


#####################################

hgnc_list <- read_csv("hgnc_list.csv")
mRNA_list <- read_csv("./2018_DEC/mRNA_list.csv")
targetList<-read_csv("targetList.csv")

mrna_finalData_selected_geneSet <- read_csv("~/Documents/Melonoma_TCGA-SKCM/mrna_finalData_selected_geneSet.csv")
mirna_significancyTestResultst <- read_csv("~/Documents/Melonoma_TCGA-SKCM/mirna_significancyTestResultst.csv")
mrna_significancyTestResultst <- read_csv("~/Documents/Melonoma_TCGA-SKCM/mrna_significancyTestResultst.csv")
mrna_finalData_selected_geneSet <- read_csv("~/Documents/Melonoma_TCGA-SKCM/mrna_finalData_selected_geneSet.csv")

# Parse version section seperated with dot.
cat(paste("\n Clone miRNA List"))
mRNA_Data <- data.frame(id = c (1:60483), rna=(1:60483), version = c(1:60483))


for (i in 1: 60483)
{
   cat("\014")
   cat(paste("..", i, ".."))
    temp<-matrix(unlist(  strsplit(mRNA_list$x[i],".", fixed=TRUE)),ncol=1, byrow=TRUE)
    mRNA_Data$id[i]<-i
    mRNA_Data$rna[i]<-temp[1,]
    mRNA_Data$version[i]<-temp[2]
    
}
View(mRNA_Data)

hgnc_list$ensembl_transcript_id<-NULL
hgnc_list$X1<-NULL
View( c)
u_Human_genes<-hgnc_list[!duplicated(hgnc_list), ]
miRTarBase_Human<-miRTarBase_MTI[miRTarBase_MTI$`Species (miRNA)`%in% 'Homo sapiens' ,]

miRNA_SET<-miRTarBase_Human[miRTarBase_Human$`Target Gene`%in% u_Human_genes$hgnc_symbol ,]
miRNA_SET$Experiments<-NULL
miRNA_SET$`Support Type`<-NULL
miRNA_SET$`References (PMID)`<-NULL
miRNA_SET<-miRNA_SET[!duplicated(miRNA_SET), ]


match.rna<-u_Human_genes[u_Human_genes$hgnc_symbol%in% miRNA_SET$`Target Gene` ,]
colnames(match.rna)<-c('mRNA','TargetGene','hgnc_id')
match.mirna<-miRNA_SET
colnames(match.mirna)<-c('mirTarbase_id','miRNA','Species_mrna','TargetGene','Target_Gene_id','Species')
df<-merge(x = match.mirna, y = match.rna, by = "TargetGene", all.x = TRUE)
# match.mirna<-match.mirna[match.mirna$miRNA %in% miRNA_list$x,]
logToFile("DatabaseOps","merged_mrna_mirna.csv",df,"D")
#write.csv(file="merged_mrna_mirna.csv",df)


# 
# ###############   PART- 2: d??zenleniyor //05.06.2019
# #################  Significannce test results collecting 
# colnames(mirna_significancyTestResultst)<-c('count','miRNA','Ftest','Ttest','Diff','P_Val','P_Val_Less','P_Val_Greater' )
# 
# 
# 
# df.mirna<-merge(x = mirna_significancyTestResultst, y = df, by = "miRNA", all.x = TRUE)
# View(df.mirna)
# 
# sample<-df.mirna[1,]
# Final.mirna_list<-sample
# View(sampleRow)
# for(i in 1:length(mirna_significancyTestResultst$miRNA))
# {
#   
#   cat("\014")
#   cat(paste("..", i, ".."))
#   part1<-mirna_significancyTestResultst[i,] 
#   phrase<-mirna_significancyTestResultst$miRNA[i]
#   sampleRow<-sample
#   sampleRow$count<-i
#   sampleRow$miRNA<-phrase
#   sampleRow$Ftest<-mirna_significancyTestResultst$Ftest[i]
#   sampleRow$Ttest<-mirna_significancyTestResultst$Ttest[i]
#   sampleRow$Diff<-mirna_significancyTestResultst$Diff[i]
#   sampleRow$P_Val<-mirna_significancyTestResultst$P_Val[i]
#   sampleRow$P_Val_Less<-mirna_significancyTestResultst$P_Val_Less[i]
#   sampleRow$P_Val_Greater<-mirna_significancyTestResultst$P_Val_Greater[i]
#   
#   search.results<- df [grep(phrase,df$miRNA ),]
#   
#   if ( length(search.results[,1])  >0)
#   {
#     for( y in 1:length(search.results[,1]) )
#     {
#       
#       part2<-search.results[y,]
#       sampleRow$TargetGene<-part2$TargetGene
#       sampleRow$mirTarbase_id<-part2$mirTarbase_id
#       sampleRow$Species_mrna<-part2$Species_mrna
#       sampleRow$Target_Gene_id<-part2$Target_Gene_id
#       sampleRow$Species<-part2$Species
#       sampleRow$mRNA<-part2$mRNA
#       sampleRow$hgnc_id<-part2$hgnc_id
#       Final.mirna_list<-rbind(Final.mirna_list,sampleRow)
#     }
#   }else
#   {
#     Final.mirna_list<-rbind(Final.mirna_list,sampleRow)
#   }
#  
# }
# 
# length(Final.mirna_list$miRNA)
# mRNA_Pval<-rep(-1,5687)
# mRNA_Diff<-rep(-1,5687)
# mRNA_Less<-rep(-1,5687)
# mRNA_Greater<-rep(-1,5687)
# Final.mirna_list<-cbind(Final.mirna_list,mRNA_Pval,mRNA_Diff,mRNA_Less,mRNA_Greater)
# View(Final.mirna_list)
# 
# for (i in 1:length(Final.mirna_list$miRNA))
# {
#   
#   cat("\014")
#   cat(paste("..", i, ".."))
#   phrase<-Final.mirna_list$mRNA[i]
#   if(!is.na(phrase))
#   {
#    
#     search.results<- mrna_significancyTestResultst [grep(phrase,mrna_significancyTestResultst$mRNA ),]
#     
#     if ( length(search.results[,1])  >0)
#     {
#       Final.mirna_list$mRNA_Pval[i]<-search.results$P_Val[1]
#       Final.mirna_list$mRNA_Diff[i]<-search.results$Diff[1]
#       Final.mirna_list$mRNA_Less[i]<-search.results$P_Val_Less[1]
#       Final.mirna_list$mRNA_Greater[i]<-search.results$P_Val_Greater[1]
#       
#     }
#   }
#   
# }
# write.csv(file="merged_mirna_significancy_results.csv",Final.mirna_list)
# 
# 
# ################# Start mRNA Aggeregation
# colnames(mrna_significancyTestResultst)<-c('count','mRNA','Ftest','Ttest','Diff','P_Val','P_Val_Less','P_Val_Greater' )
# 
# df.mrna<-merge(x = mrna_significancyTestResultst, y = df, by = "mRNA", all.x = TRUE)
# View(df.mrna)
# 
# 
# mrna_significancyTestResultst_clone<-mrna_significancyTestResultst[mrna_significancyTestResultst$Diff %in% 1,]
# sample<-df.mrna[1,]
# Final.mrna_list<-sample
# for(i in 1:length(mrna_significancyTestResultst_clone$mRNA))
# {
#   cat("\014")
#   cat(paste("..", i, ".."))
#   part1<-mrna_significancyTestResultst_clone[i,] 
#   phrase<- as.vector(unlist( strsplit(mrna_significancyTestResultst_clone$mRNA[i],".",fixed=TRUE)))[1]
#   sampleRow<-sample
#   sampleRow$count<-i
#   sampleRow$mRNA<-mrna_significancyTestResultst_clone$mRNA[i]
#   sampleRow$Ftest<-mrna_significancyTestResultst_clone$Ftest[i]
#   sampleRow$Ttest<-mrna_significancyTestResultst_clone$Ttest[i]
#   sampleRow$Diff<-mrna_significancyTestResultst_clone$Diff[i]
#   sampleRow$P_Val<-mrna_significancyTestResultst_clone$P_Val[i]
#   sampleRow$P_Val_Less<-mrna_significancyTestResultst_clone$P_Val_Less[i]
#   sampleRow$P_Val_Greater<-mrna_significancyTestResultst_clone$P_Val_Greater[i]
#   
#   search.results<- df [grep(phrase,df$mRNA ),]
#   cat(phrase)
#   if ( length(search.results[,1])  >0)
#   {
#     for( y in 1:length(search.results[,1]) )
#     {
#       cat(paste("..", i, ".",y))
#       part2<-search.results[y,]
#       sampleRow$TargetGene<-part2$TargetGene
#       sampleRow$mirTarbase_id<-part2$mirTarbase_id
#       sampleRow$Species_mrna<-part2$Species_mrna
#       sampleRow$Target_Gene_id<-part2$Target_Gene_id
#       sampleRow$Species<-part2$Species
#       sampleRow$miRNA<-part2$miRNA
#       sampleRow$hgnc_id<-part2$hgnc_id
#       Final.mrna_list<-rbind(Final.mrna_list,sampleRow)
#     }
#   }else
#   {
#     Final.mrna_list<-rbind(Final.mrna_list,sampleRow)
#   }
#   
# }
# 
# length(Final.mrna_list$mRNA)
# miRNA_Pval<-rep(-1,9665)
# miRNA_Diff<-rep(-1,9665)
# miRNA_Pval_Less<-rep(-1,9665)
# miRNA_Pval_Greater<-rep(-1,9665)
# Final.mrna_list<-cbind(Final.mrna_list,miRNA_Pval,miRNA_Diff,miRNA_Pval_Less,miRNA_Pval_Greater)
# View(Final.mrna_list)
# Final.mrna_list<-merged_mrna_significancy_results
# 
# for (i in 1:length(Final.mrna_list$mRNA))
# {
# 
#   cat("\014")
#   cat(paste("..", i, ".."))
#   phrase<-Final.mrna_list$miRNA[i]
#   if(!is.na(phrase))
#   {
#     search.results<- mirna_significancyTestResultst [grep(sapply(phrase,tolower),mirna_significancyTestResultst$miRNA ),]
#     
#     if ( length(search.results[,1])  >0)
#     {
#       Final.mrna_list$miRNA_Pval[i]<-search.results$P_Val[1]
#       Final.mrna_list$miRNA_Diff[i]<-search.results$Diff[1]
#       Final.mrna_list$miRNA_Pval_Less[i]<-search.results$P_Val_Less[1]
#       Final.mrna_list$miRNA_Pval_Greater[i]<-search.results$P_Val_Greater[1]
#       
#     }
#   }
#   
# }
# 
# write.csv(file="merged_mrna_significancy_results.csv",Final.mrna_list)
# 
# ############################################################################
# 
# Sig.mRNA <- Final.mrna_list[Final.mrna_list$Diff %in% 1, ]
# Sig.mRNA_sigMirna <- Sig.mRNA[Sig.mRNA$miRNA_Diff %in% 1,]
# 
# Sig.miRNA <- Final.mirna_list[Final.mirna_list$P_Val < 0.05, ]
# Sig.miRNA_SigMrna <- Sig.miRNA[Sig.miRNA$mRNA_Diff %in% 1,]
# 
# 
# Sig.mRNA_sigMirna$X1<-NULL
# Sig.mRNA_sigMirna$count<-NULL
# Sig.mRNA_sigMirna$Ftest<-NULL
# Sig.mRNA_sigMirna$Ttest<-NULL
# Sig.mRNA_sigMirna$TargetGene<-NULL
# Sig.mRNA_sigMirna$TargetGene<-NULL
# Sig.mRNA_sigMirna$mirTarbase_id<-NULL
# Sig.mRNA_sigMirna$Species<-NULL
# Sig.mRNA_sigMirna$hgnc_id<-NULL
# Sig.mRNA_sigMirna$Species_mrna<-NULL
# Sig.mRNA_sigMirna$Target_Gene_id<-NULL
# 
# Sig.miRNA_SigMrna$X1<-NULL
# Sig.miRNA_SigMrna$count<-NULL
# Sig.miRNA_SigMrna$Ftest<-NULL
# Sig.miRNA_SigMrna$Ttest<-NULL
# Sig.miRNA_SigMrna$TargetGene<-NULL
# Sig.miRNA_SigMrna$TargetGene<-NULL
# Sig.miRNA_SigMrna$mirTarbase_id<-NULL
# Sig.miRNA_SigMrna$Species<-NULL
# Sig.miRNA_SigMrna$hgnc_id<-NULL
# Sig.miRNA_SigMrna$Species_mrna<-NULL
# Sig.miRNA_SigMrna$Target_Gene_id<-NULL
# 
# colnames(Sig.miRNA_SigMrna)<-c("prob","diff","P_Val","P_Val_less","P_Val_Greater","target","target_P_Val","target_diff","target_P_Val_less","target_P_Val_Greater")
# colnames(Sig.mRNA_sigMirna)<-c("prob","diff","P_Val","P_Val_less","P_Val_Greater","target","target_P_Val","target_diff","target_P_Val_less","target_P_Val_Greater")
# 
# 
# 
# for ( i in 1:34)
# {
#   Sig.miRNA_SigMrna$miRNA[i] <- 0
#   if (Sig.miRNA_SigMrna$P_Val_less[i]>Sig.miRNA_SigMrna$P_Val_Greater[i])
#   {
#     Sig.miRNA_SigMrna$miRNA[i] <- 1
#   }else if (Sig.miRNA_SigMrna$P_Val_less[i]<Sig.miRNA_SigMrna$P_Val_Greater[i])
#     
#   {
#     Sig.miRNA_SigMrna$miRNA[i] <- -1
#   }
#   
#   if (Sig.miRNA_SigMrna$target_P_Val_less[i]>Sig.miRNA_SigMrna$target_P_Val_Greater[i])
#   {
#     Sig.miRNA_SigMrna$mRNA[i] <- 1 
#   }else if (Sig.miRNA_SigMrna$target_P_Val_less[i]<Sig.miRNA_SigMrna$target_P_Val_Greater[i])
#     
#   {
#     Sig.miRNA_SigMrna$mRNA[i] <- -1
#   }
#   
# }
# 
# 
# 
# for ( i in 1:122)
# {
#   
#   if (Sig.mRNA_sigMirna$P_Val_less[i]>Sig.mRNA_sigMirna$P_Val_Greater[i])
#   {
#     Sig.mRNA_sigMirna$mRNA[i] <- 1 
#   }else if (Sig.mRNA_sigMirna$P_Val_less[i]<Sig.mRNA_sigMirna$P_Val_Greater[i])
#     
#   {
#     Sig.mRNA_sigMirna$mRNA[i] <- -1 
#   }
#   
#   if (Sig.mRNA_sigMirna$target_P_Val_less[i]>Sig.mRNA_sigMirna$target_P_Val_Greater[i])
#   {
#     Sig.mRNA_sigMirna$miRNA[i] <- 1 
#   }else if (Sig.mRNA_sigMirna$target_P_Val_less[i]<Sig.mRNA_sigMirna$target_P_Val_Greater[i])
#     
#   {
#     Sig.mRNA_sigMirna$miRNA[i] <- -1
#   }
# }
# 
# 
# selected <- rbind(Sig.mRNA_sigMirna, Sig.miRNA_SigMrna)
# write.csv(file="final_selected_set.csv",selected)





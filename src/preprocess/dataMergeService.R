margeData <-function()
{
  
  path <- paste(global.workingDir,"/data/raw/",sep = "")
  miRNA_data_final <- read_csv(paste(path,global.datafile.miRna.raw,sep=""))
  htseq_data_final <- read_csv(paste(path,global.datafile.mRna.raw,sep=""))
  
  nframe_rna<-subset(htseq_data_final, select= -c(1:20,22:24))
  nframe_mirna<-miRNA_data_final
  
  df<-merge(x = nframe_mirna, y = nframe_rna, by = "case_id_1", all = TRUE)
  df.notNull<-removeNullClass(df)
  df.notNullNotZero<- removeZeroVariables(df.notNull,25)
  
  colnames(df.notNullNotZero)[11]<-"Class"
  finalDataSet<-subset(df.notNullNotZero, select = -c(2:10,12:24))
  finalDataSet<-na.omit(finalDataSet)
  finalDataSet

}
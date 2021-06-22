

#for the latest version install from github

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
source("https://bioconductor.org/biocLite.R")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")


library("TCGAbiolinks")
library("SummarizedExperiment")


# Download and format TCGA OV microarray data
#try and download the microarray expressions data.
query_microarray <- GDCquery(project = "TCGA-OV", 
                             data.category = "Gene expression",
                             data.type = "Gene expression quantification",
                             platform = "HT_HG-U133A",
                             access = "open",
                             legacy = TRUE)
GDCdownload(query_microarray )
OVMicroarray <- GDCprepare(query_microarray )

microarray <- assay(OVMicroarray)

#remove the duplicate genes and make gene names the matrix row names
microarray <- microarray[which(!duplicated(rownames(microarray))),]

#compute the 12 character barcode for each patients
microarrayPatients <- cbind(colnames(microarray), gsub('\\.','-',substring(colnames(microarray),1,12)))

#only include patients that were included in Verhaak dataset
# microarray <- microarray[,which(microarrayPatients[,2] %in% classDefinitions_verhaak[which(!is.na(classDefinitions_verhaak$SUBTYPE)),"ID"])]


microarrayPatients <- merge(microarrayPatients,classDefinitions_verhaak[,c("ID","SUBTYPE")],by.x = 2, by.y =1)
colnames(microarrayPatients) <- c( "barcode","patient","SUBTYPE")

#only include patients that have microarray data for them
microarrayPatients <- microarrayPatients[which(microarrayPatients$patient %in% colnames(microarray)),]
microarrayPatients <- microarrayPatients[order(microarrayPatients$SUBTYPE),]

#convert the barcodes so that they will be compatible with colnames (R doesn't like "-" in column names)
microarrayPatients$patient <- gsub('-','\\.',microarrayPatients$patient)
colnames(microarray) <- gsub('-','\\.',colnames(microarray))
microarray <- microarray[,colnames(microarray)[order(match(colnames(microarray),microarrayPatients$patient))]]
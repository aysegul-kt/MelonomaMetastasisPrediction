p_vals <- c(0.001)
base_item <- c(2,3)
definition_items <-  c(6) # rbf eklenecek...
definition_list <- c( "Type_2_3_6")  # rbf eklenecek...

typeList <- c( base_item,  definition_items[t])
selectedVar <- getSignificantVariables(flaggedSigResult, p_vals[1],typeList)

NetworkData <- flaggedSigResult[flaggedSigResult$relation_type %in% c(2,3,6),]

logToFile("Network","group_2_3_6.csv",NetworkData , "D")

mRNA_unique_list <- unique(NetworkData$targetName)

NetworkData$targetGene <- rep(c( NA),length(NetworkData[,1]))
NetworkData$targetGeneId <- rep(c( NA),length(NetworkData[,1]))
NetworkData$hgnc_id <- rep(c( NA),length(NetworkData[,1]))
NetworkData$ensembl_gene_id <- rep(c( NA),length(NetworkData[,1]))


i<-1

for (i in 1:length(mRNA_unique_list) )
{
  cat ( paste ( "\n ", i ))
   parsed <- unlist(strsplit( as.character(  mRNA_unique_list[i]),".", fixed=TRUE))[1]
   selected_def <- mrna_mirna_mapping [mrna_mirna_mapping$mRNA %in% parsed ,]  # One row is enhough to get definition since all are equal
   rcount <- nrow ( NetworkData[NetworkData$targetName %in% mRNA_unique_list[i],])
    if ( nrow(selected_def)>0)
    {
     NetworkData$targetGene[NetworkData$targetName %in% mRNA_unique_list[i] ] <- rep(c( selected_def$TargetGene[1]),rcount)
     NetworkData$targetGeneId[NetworkData$targetName %in% mRNA_unique_list[i] ] <- rep(c( selected_def$Target_Gene_id[1]),rcount)
     NetworkData$hgnc_id[NetworkData$targetName %in% mRNA_unique_list[i] ] <- rep(c( selected_def$hgnc_id[1]),rcount)
     NetworkData$ensembl_gene_id[NetworkData$targetName %in% mRNA_unique_list[i] ] <- rep(c( selected_def$mRNA[1]),rcount)
    }
}




miRNAPart <- NetworkData[NetworkData$ProbName %in% selectedVar , ]
mRNAPart <- NetworkData[NetworkData$targetName %in% selectedVar , ]
combined_network  <- distinct(rbind(mRNAPart, miRNAPart))

library(dplyr)

View(combined_network)

logToFile("Network","mRNA_gene_list.csv",unique(combined_network$targetGene) , "D")

logToFile("Network","selectedVariables.csv",selectedVar , "D")

logToFile("Network","combined_network.csv",combined_network , "D")





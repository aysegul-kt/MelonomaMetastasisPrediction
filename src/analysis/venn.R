library(RColorBrewer)
library(VennDiagram)


myCol <- brewer.pal(2, "Pastel1")

set1 <- unlist(lapply(as.vector( unique(ContributionsFrame.Step2$attribute)),  function(x) {str_replace_all(x,"_x_","-")}))
set2 <- unlist(lapply(as.vector( unique(ContributionsFrame.Step3$attribute)),  function(x) {str_replace_all(x,"_x_","-")}))

lvariables <- rbind(pool$mRNAFeatures,pool$methylationFeatures)
setorder(lvariables,P_Val) 

lvariablesList <- as.vector(lvariables[1:2000,]$probName)  
set3 <- as.vector( unique(lvariablesList))

venn.diagram(
  x = list(set1, set2),
  category.names = c("Set 1" , "Set 2 " ),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)


commonSet <- unique(set2[(set2 %in% set1)])


diffSet <- setdiff(set2, set1) 


View(diffSet)


temp  <-  getReleatedMarkers(set2, NULL)
temp2 <-  getReleatedMarkers(diffSet, NULL)
temp3  <-  getReleatedMarkers(set1, NULL)
temp_set3 <- getReleatedMarkers(set3, NULL)

View(temp)
View(temp2)
View(temp3)
View(temp_set3)

logToFile("Venn",paste("SelectedMarkersCrossTable_all.csv"), temp ,"L")
logToFile("Venn",paste("SelectedMarkersCrossTable_diffSet.csv"), temp2 ,"L")
logToFile("Venn",paste("SelectedMarkersCrossTable_all_step2.csv"), temp3 ,"L")
logToFile("Venn",paste("SelectedMarkersCrossTable_all_step4.csv"), temp_set3 ,"L")




##############################


temp_list <- unlist(strsplit(davidoutput_step3$`EnrichmentMap::Genes`[1],"[|]")[[1]])

pathway_genes <- data.frame(ensembl_gene_id  <-   temp_list, pathways <- rep(davidoutput_step3$name[1],length(temp_list)  ),stringsAsFactors=FALSE)

names(pathway_genes) <- c("ensembl_gene_id","pathways")

for (i in c(2:14)) {
  
  temp_genes <- unlist(strsplit(davidoutput_step3$`EnrichmentMap::Genes`[i],"[|]")[[1]])
  
  for (gene in temp_genes) {
   cat( gene)
   dummyrow <-  pathway_genes[pathway_genes$ensembl_gene_id %in% gene,]
  
   if (nrow( dummyrow) >0){
     cat("append")
     pathway_genes[pathway_genes$ensembl_gene_id %in% gene,]$pathways <-  paste(dummyrow$pathways,"|",davidoutput_step3$name[i],sep = "")
   }else{
     cat("new")
     pathway_genes <-  rbind(pathway_genes, c(gene,davidoutput_step3$name[i] ))
   }
  }
}

View(pathway_genes)
result <- merge(Venn_SelectedMarkersCrossTable_all_cytoscapeInput_2,pathway_genes , by = "ensembl_gene_id",all= TRUE, x.all=TRUE)
result <- merge(result,pathway_genes , by = "ensembl_gene_id",all= TRUE, x.all=TRUE)


logToFile("Venn",paste("SelectedMarkersCrossTable_all_cytoscapeInput_pathwaymerged.csv"), result ,"L")


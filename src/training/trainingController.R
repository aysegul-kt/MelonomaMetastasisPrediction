getSignificantVariables  <- function(result, p_val, selectedset){
  
  selectedVars <- c(NA)
 
   for (i in 1:nrow(result) )
  {
    #cat("\014")
    #cat(i)
    if (result$relation_type[i] %in% selectedset)  
    {
      if (!is.na(result$ProbName[i] ) ){
        if (result$P_Val[i] < p_val) { selectedVars <- rbind (selectedVars,as.character(result$ProbName[i])) }
      }
 
      if (!is.na(result$targetName[i])) {
        if(result$targetP_Val[i] < p_val)  { selectedVars <- rbind(selectedVars,as.character(result$targetName[i])) }
      } 
    }
    
  }
  
  s<-(unique(na.omit(selectedVars)))
  #cat("\014")
  cat(length(s))
  
  
  # mySet<-result[result$relation_type %in% selectedset, ]
  s

}

FilterSignificantVariables  <- function(FilterInput, includeMety){


miRNA <- FilterInput[FilterInput$miRNA_Diff001 %in% "1",]

mRNA <- FilterInput[FilterInput$mRNA_Diff001 %in% "1",]

V1 <- as.vector( miRNA$miRNA_probName)
V2 <- as.vector( mRNA$mRNA_probName)
V3 <- as.vector( mRNA$mRNA_probName)


if (includeMety) {
  methy <- FilterInput[FilterInput$methy_Diff01 %in% "1",] # p val = 0.001 stored on diff column.... for methylation...
  V3 <- as.vector( methy$methy_probName)
 
}

selection <- c( V1,V2,V3)

 selection<-(unique(na.omit(selection)))
 #cat("\014")
  cat(length(selection))

  selection
  
}
dataNormalization <- function(startAt, tempDataSet, stopAt)
{
startindex<-startAt
if(is.null( stopAt)) {
  stopAt <- ncol(tempDataSet)
}

for (i in startindex:ncol(tempDataSet)){
  cat(paste("\n", i , "of", ncol(tempDataSet)))
  tempDataSet[i]<- normalize(tempDataSet[i], method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
  cat("\014")
}
tempDataSet
}

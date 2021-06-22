logToFile <- function (ParentClazz, fileName, data, logtype)
{
  if(logtype == "L" ){write.csv(file=paste(global.workingDir,"/output/logs/",ParentClazz,"_", fileName, sep = ""),data)}
  else  if (logtype == "F" ){ write.csv(file=paste(global.workingDir,"/output/figures/",ParentClazz,"_", fileName, sep = ""),data)}
  else if  (logtype == "D" ) { write.csv(file=paste(global.workingDir,"/output/IntermediateData/",ParentClazz,"_", fileName, sep = ""),data)}
  else {write.csv(file=paste(global.workingDir,"/output/logs/",ParentClazz,"_", fileName, sep = ""),data)}
  
}

ReadFromLogFile <- function (ParentClazz, fileName, logtype)
{
  if(logtype == "L" ){
    data <-  read.csv(file=paste(global.workingDir,"/output/logs/",ParentClazz,"_", fileName, sep = ""))
    }
  else  if (logtype == "F" ){
    data <- read.csv(file=paste(global.workingDir,"/output/figures/",ParentClazz,"_", fileName, sep = ""))
    }
  else if  (logtype == "D" ) { 
    data <- read.csv(file=paste(global.workingDir,"/output/IntermediateData/",ParentClazz,"_", fileName, sep = ""))
    }
  else {
    data <- read.csv(file=paste(global.workingDir,"/output/logs/",ParentClazz,"_", fileName, sep = ""))
    }
  
   
}
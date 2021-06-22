###############################################
# ADI              :  Base An. Script
# YAZAR            : AYSEGUL KUTLAY
# ACIKLAMA         :
# GUNCELLEME TAR?H?: 06.06.2019
# VERSIYON         : 02
################################################
TrimMiRNA <- function(trimInputdata)
{
  trimInputdata <- joind_map
  trimInputdata$refName  <- rep("", nrow(trimInputdata))
  tempMirnaList <- unique(trimInputdata$miRNA)
  
  for (i in 1:length(tempMirnaList))
  {
    mirnaname <- str_to_lower(tempMirnaList[i])
    condition <- str_count(mirnaname, '-')
    
    while (condition>2) {
      mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 1)
      condition <- str_count(mirnaname, '-')
    }
   

    # x$Age[x$Name %in% "John"] <- rep(c(10),length(x$Age[x$Name %in% "John"] )) 
    trimInputdata$refName[trimInputdata$miRNA %in% tempMirnaList[i]] <- rep(c(mirnaname),length(trimInputdata$refName[trimInputdata$miRNA %in% tempMirnaList[i]]))
    
    cat(paste(i, ": ", tempMirnaList[i], " to ", mirnaname, "\n"))
  }
  
  trimInputdata
}

getRegulationDetailsResults <- function(mydata)
{
  mydata$direction <- rep(NA, nrow(mydata))
  mydata$Diff001 <- rep(c(1), nrow(mydata))
  mydata$Diff005 <- rep(c(1), nrow(mydata))
  mydata$dtype <- rep(NA, nrow(mydata))
  mydata$refName <- rep("", nrow(mydata))
  
  for (i in 1:nrow(mydata))
  {
    if (mydata$P_Val_Less[i] > mydata$P_Val_Greater[i])
    {
      mydata$direction[i] <- "H"
    } else
    {
      mydata$direction[i] <- "L"
    }
    
    if (mydata$P_Val[i] > 0.001) {
      mydata$Diff001[i] <- 0
    }
    if (mydata$P_Val[i] > 0.005) {
      mydata$Diff005[i] <- 0
    }
    
    if (i < 1608) {
      mydata$dtype[i] <- "miRNA"
      mirnaname <- str_to_lower(mydata$probName[i])
      if (str_count(mirnaname, '-') == 3) {
        mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 2)
        if (str_count(mirnaname, '-') == 3) {
          mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 1)
          if (str_count(mirnaname, '-') == 3) {
            mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 1)
            
          }
        }
      }
      mydata$refName[i] <- mirnaname
    }
    else {
      mydata$dtype[i] <- "mRNA"
      phrase <- mydata$probName[i]
      phrase <- str_sub(phrase, 0, str_locate(phrase, "\\.")[1] - 1)
      mydata$refName[i] <- phrase
    }
    
    
    
    
    
  }
  
  mydata
  
}

getGeneMethylationDetailedResults <- function(mydata){
  mydata$direction <- rep(NA, nrow(mydata))
  mydata$Diff001 <- rep(c(1), nrow(mydata))
  mydata$Diff005 <- rep(c(1), nrow(mydata))
  mydata$dtype <- rep("Methylation", nrow(mydata))
  mydata$hgnc_symbol <- mydata$ProbName
  
  for (i in 1:nrow(mydata))
  {
    if (mydata$P_Val_Less[i] > mydata$P_Val_Greater[i])
    {
      mydata$direction[i] <- "H"
    } else
    {
      mydata$direction[i] <- "L"
    }
    
    if (mydata$P_Val[i] > 0.0001) {
      mydata$Diff001[i] <- 0
    }
    
    if (mydata$P_Val[i] > 0.0005) {
      mydata$Diff005[i] <- 0
    }
    
 
  }
  
  mydata
}
getPattern <- function(miRNAPatern, mRNAPatern, methylationPatern, grp_code)
{
  # uses
  # ttest.results_detailed
  # joind_map
  # miRNAPatern <-  "NA"
  # mRNAPatern <- "H"
  # grp_code <-"6"
  
 # NOTE: null target of mRNA will be added 28.02.2020
  
  subset.miRNA.all <-
    ttest.results_detailed[ttest.results_detailed$dtype %in% "miRNA",]
  
  if (miRNAPatern != "NA") {
    subset.miRNA.pattern <-
      subset.miRNA.all[subset.miRNA.all$Diff001  %in% "1", ]
    subset.miRNA.pattern <-
      subset.miRNA.pattern[subset.miRNA.pattern$direction  %in% miRNAPatern, ]
  } else if (miRNAPatern == "NA"  ){
    subset.miRNA.pattern <-
      subset.miRNA.all[subset.miRNA.all$Diff001  %in% "0", ]
  }else {
    subset.miRNA.pattern <-
      subset.miRNA.all[subset.miRNA.all$Diff001  %in% "nan", ]
  }
 
  #View(subset.miRNA.pattern)
  
  
  subset.mRNA.all <-
    ttest.results_detailed[ttest.results_detailed$dtype %in% "mRNA",]
  
  if (mRNAPatern != "NA"  ) {
    subset.mRNA.pattern <-
      subset.mRNA.all[subset.mRNA.all$Diff001  %in% "1", ]
    subset.mRNA.pattern <-
      subset.mRNA.pattern[subset.mRNA.pattern$direction  %in% mRNAPatern, ]
  } else if (mRNAPatern == "NA"  ){
    subset.mRNA.pattern <- subset.mRNA.all[subset.mRNA.all$Diff001  %in% "0", ]
  }else {
    subset.mRNA.pattern <- subset.mRNA.all[subset.mRNA.all$Diff001  %in% "nan", ]
  }
  
  if (methylationPatern != "Nan") {
    subset.methylation.pattern <- GeneMethylationDetailedResults[GeneMethylationDetailedResults$Diff001  %in% "1", ]
    subset.methylation.pattern <- subset.methylation.pattern[subset.methylation.pattern$direction  %in% methylationPatern, ]
  } else{
    subset.methylation.pattern <-  GeneMethylationDetailedResults[GeneMethylationDetailedResults$Diff001  %in% c("nan"), ]
  }
  #View(subset.mRNA.pattern)
  #View(joind_map)
 
  miRNA.Set <- subset.miRNA.pattern
  #View(miRNA.Set)
  
  subset.joinedmap <-joind_map[joind_map$refName %in% subset.miRNA.pattern$refName,]
 
  #View(subset.joinedmap)
  myList <- unique(subset.joinedmap$ensembl_gene_id)
  
  mRNA.Set <-subset.mRNA.pattern[subset.mRNA.pattern$refName %in%  myList, ]
  #View(mRNA.Set)
  
  
  
  colnames(miRNA.Set) <-
    c(
      "miRNA_probName",
      "miRNA_Ftest",
      "miRNA_Ttest",
      "miRNA_Diff01",
      "miRNA_P_Val",
      "miRNA_P_Val_Less",
      "miRNA_P_Val_Greater" ,
      "miRNA_direction",
      "miRNA_Diff001",
      "miRNA_Diff005",
      "miRNA_dtype",
      "refName"
    )
  colnames(mRNA.Set) <-
    c(
      "mRNA_probName",
      "mRNA_Ftest",
      "mRNA_Ttest",
      "mRNA_Diff01",
      "mRNA_P_Val",
      "mRNA_P_Val_Less",
      "mRNA_P_Val_Greater" ,
      "mRNA_direction",
      "mRNA_Diff001",
      "mRNA_Diff005",
      "mRNA_dtype",
      "ensembl_gene_id"
    )
  
  
  
    #View(subset.joinedmap)
    myList <- unique(subset.joinedmap$`Target Gene`)
    
    # `Target Gene`
    methylation.Set <-   subset.methylation.pattern[subset.methylation.pattern$`Target Gene` %in%  myList, ]
    
    colnames(methylation.Set) <-
      c(
        "methy_probName",
        "methy_Ftest",
        "methy_Ttest",
        "methy_Diff01",
        "methy_P_Val",
        "methy_P_Val_Less",
        "methy_P_Val_Greater" ,
        "methy_direction",
        "methy_Diff001",
        "methy_Diff005",
        "methy_dtype",
        "hgnc_symbol"
      )
  
  
  
  
  
  result <- merge(subset.joinedmap, miRNA.Set, by = "refName")

  result <- merge(result, mRNA.Set, by = "ensembl_gene_id")
  
  result <- merge(x=result, y=methylation.Set, by = "hgnc_symbol",all= TRUE, x.all=TRUE)
  
    
  #View(result)
  # colnames(result[    c(1, 5, 8, 9, 10, 11, 12, 13, 14, 18, 19, 29, 30,40,41)])
  
  result <- result[, - c(1, 5, 8, 9, 10, 11, 12, 13, 14, 18, 19, 29, 30,40,41)]
  result <- result[!duplicated(result), ]
  
  result$grp_code <- rep(c(grp_code), nrow(result))
  
  result
  
}

getReleatedMarkers <- function(markerList,grp_code)
{
  # uses
  # ttest.results_detailed
  # joind_map
  # miRNAPatern <-  "NA"
  # mRNAPatern <- "H"
  # grp_code <-"6"
  
  # NOTE: null target of mRNA will be added 28.02.2020
  
  subset.miRNA.all <-
    ttest.results_detailed[ttest.results_detailed$dtype %in% "miRNA",]
  

  miRNA.Set <-
      subset.miRNA.all[subset.miRNA.all$probName  %in%  markerList, ]
    
   subset.mRNA.all <-
    ttest.results_detailed[ttest.results_detailed$dtype %in% "mRNA",]
  
  mRNA.Set <-
      subset.mRNA.all[subset.mRNA.all$probName  %in% markerList, ]
  
  methylation.Set <- GeneMethylationDetailedResults[GeneMethylationDetailedResults$ProbName  %in% markerList, ]
  
  subset.joinedmap <-rbind(joind_map[joind_map$refName %in% miRNA.Set$refName,],
                           joind_map[joind_map$ensembl_gene_id %in% miRNA.Set$refName,],
                           joind_map[joind_map$hgnc_symbol %in% methylation.Set$ProbName,])
                           
  
  
  
  colnames(miRNA.Set) <-
    c(
      "miRNA_probName",
      "miRNA_Ftest",
      "miRNA_Ttest",
      "miRNA_Diff01",
      "miRNA_P_Val",
      "miRNA_P_Val_Less",
      "miRNA_P_Val_Greater" ,
      "miRNA_direction",
      "miRNA_Diff001",
      "miRNA_Diff005",
      "miRNA_dtype",
      "refName"
    )
  colnames(mRNA.Set) <-
    c(
      "mRNA_probName",
      "mRNA_Ftest",
      "mRNA_Ttest",
      "mRNA_Diff01",
      "mRNA_P_Val",
      "mRNA_P_Val_Less",
      "mRNA_P_Val_Greater" ,
      "mRNA_direction",
      "mRNA_Diff001",
      "mRNA_Diff005",
      "mRNA_dtype",
      "ensembl_gene_id"
    )
  
  
  
  
  colnames(methylation.Set) <-
    c(
      "methy_probName",
      "methy_Ftest",
      "methy_Ttest",
      "methy_Diff01",
      "methy_P_Val",
      "methy_P_Val_Less",
      "methy_P_Val_Greater" ,
      "methy_direction",
      "methy_Diff001",
      "methy_Diff005",
      "methy_dtype",
      "hgnc_symbol"
    )
  
  
  
  
  
  result <- merge(subset.joinedmap,miRNA.Set , by = "refName")
  
  result <- merge(result, mRNA.Set, by = "ensembl_gene_id",all= TRUE, x.all=TRUE)
  
  result <- merge(x=result, y=methylation.Set, by = "hgnc_symbol",all= TRUE, x.all=TRUE)
  
  
  #View(result)
  # colnames(result[    c(1, 5, 8, 9, 10, 11, 12, 13, 14, 18, 19, 29, 30,40,41)])
  
  result <- result[, - c(1, 5, 8, 9, 10, 11, 12, 13, 14, 18, 19, 29, 30,40,41)]
  result <- result[!duplicated(result), ]
  
  result$grp_code <- rep(c(grp_code), nrow(result))
  
  result
  
}

getCrossTableForSignificancyResults <- function()
{
  path <-
    paste(global.workingDir,
          "/data/raw/",
          global.mrna_mirna_mappingFile,
          sep = "")
  
  mrna_mirna_mapping <- read_csv(path)
  significancyTestResults <- ttest.results
  
  
  result.part0 <-
    data.frame(
      ProbName = character(0),
      Ftest = character(0),
      Ttest = character(0),
      Diff = numeric(0),
      P_Val = numeric(0),
      P_Val_Less = numeric(0),
      P_Val_Greater = numeric(0),
      targetName = character(0),
      targetFtest = character(0),
      targetTtest = character(0),
      targetDiff = numeric(0),
      targetP_Val = numeric(0),
      targetP_Val_Less = numeric(0),
      targetP_Val_Greater = numeric(0),
      relation_type = character(0)
    )
  
  sample <- result.part0[1, ]
  
  for (i in 1:1607)
  {
    cat("\014")
    cat(paste("..", i, ".."))
    sample$ProbName <- significancyTestResults$probName[i]
    sample$Ftest <- significancyTestResults$Ftest[i]
    sample$Ttest <- significancyTestResults$Ttest[i]
    sample$Diff <- significancyTestResults$Diff[i]
    sample$P_Val <- significancyTestResults$P_Val[i]
    sample$P_Val_Less <- significancyTestResults$P_Val_Less[i]
    sample$P_Val_Greater <- significancyTestResults$P_Val_Greater[i]
    result.part0 <- rbind(result.part0, sample)
  }
  
  result.part1 <-
    runmRNA(1608, 2000, significancyTestResults, mrna_mirna_mapping)
  result.part2 <-
    runmRNA(2001, 5000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part2.csv",
            result.part2,
            "L")
  
  result.part3 <-
    runmRNA(5001, 10000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part3.csv",
            result.part3,
            "L")
  
  result.part4 <-
    runmRNA(10001, 15000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part4.csv",
            result.part4,
            "L")
  
  result.part5 <-
    runmRNA(15001, 20000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part5.csv",
            result.part5,
            "L")
  
  result.part6 <-
    runmRNA(20001, 25000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part6.csv",
            result.part6,
            "L")
  
  result.part7 <-
    runmRNA(25001, 30000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part7.csv",
            result.part7,
            "L")
  
  result.part8 <-
    runmRNA(30001, 35000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part8.csv",
            result.part8,
            "L")
  
  result.part9 <-
    runmRNA(35001, 40000, significancyTestResults, mrna_mirna_mapping)
  logToFile("SignificancyDataCol",
            "result.part9.csv",
            result.part9,
            "L")
  
  result.part10 <-
    runmRNA(
      40001,
      nrow(significancyTestResults),
      significancyTestResults,
      mrna_mirna_mapping
    )
  
  
  result.part3 <-
    read_csv("output/logs/SignificancyDataCol_result.part3.csv",
             col_types = cols(X1 = col_skip()))
  result.part4 <-
    read_csv("output/logs/SignificancyDataCol_result.part4.csv",
             col_types = cols(X1 = col_skip()))
  result.part5 <-
    read_csv("output/logs/SignificancyDataCol_result.part5.csv",
             col_types = cols(X1 = col_skip()))
  result.part6 <-
    read_csv("output/logs/SignificancyDataCol_result.part6.csv",
             col_types = cols(X1 = col_skip()))
  result.part7 <-
    read_csv("output/logs/SignificancyDataCol_result.part7.csv",
             col_types = cols(X1 = col_skip()))
  result.part8 <-
    read_csv("output/logs/SignificancyDataCol_result.part8.csv",
             col_types = cols(X1 = col_skip()))
  result.part9 <-
    read_csv("output/logs/SignificancyDataCol_result.part9.csv",
             col_types = cols(X1 = col_skip()))
  
  
  
  r <-
    rbind(
      result.part0,
      result.part1,
      result.part2,
      result.part3,
      result.part4,
      result.part5,
      result.part6,
      result.part7,
      result.part8,
      result.part9,
      result.part10
    )
  r
}


runmRNA <-
  function(startAt,
           StopAt,
           significancyTestResults,
           mrna_mirna_mapping)
  {
    result.part <-
      data.frame(
        ProbName = character(0),
        Ftest = character(0),
        Ttest = character(0),
        Diff = numeric(0),
        P_Val = numeric(0),
        P_Val_Less = numeric(0),
        P_Val_Greater = numeric(0),
        targetName = character(0),
        targetFtest = character(0),
        targetTtest = character(0),
        targetDiff = numeric(0),
        targetP_Val = numeric(0),
        targetP_Val_Less = numeric(0),
        targetP_Val_Greater = numeric(0),
        relation_type = character(0)
      )
    
    sample <- result.part[1, ]
    subset <- significancyTestResults
    # subset<-significancyTestResults[significancyTestResults$Diff==1, ]
    
    
    for (i in startAt:StopAt)
      #nrol(significancyTestResults))
    {
      cat("\014")
      cat(paste("..", i, ".."))
      
      # create mrna Part
      phrase <- significancyTestResults$probName[i]
      
      phrase <- str_sub(phrase, 0, str_locate(phrase, "\\.")[1] - 1)
      
      sample$targetName <- significancyTestResults$probName[i]
      sample$targetFtest <- significancyTestResults$Ftest[i]
      sample$targetTtest <- significancyTestResults$Ttest[i]
      sample$targetDiff <- significancyTestResults$Diff[i]
      sample$targetP_Val <- significancyTestResults$P_Val[i]
      sample$targetP_Val_Less <- significancyTestResults$P_Val_Less[i]
      sample$targetP_Val_Greater <-
        significancyTestResults$P_Val_Greater[i]
      
      result.part <- rbind(result.part, sample)
      
      #search target mRNA of given mirna
      search.results <-
        na.omit(mrna_mirna_mapping [grep(phrase, mrna_mirna_mapping$mRNA), ])
      
      if (nrow(search.results)  > 0)
      {
        for (y in 1:nrow(search.results))
        {
          cat("\014")
          cat(paste("..", i, "..", nrow(search.results) , "/", y))
          # for each  mrna search if there are any significacy result
          mirnaname <- str_to_lower(search.results$miRNA[y])
          if (str_count(mirnaname, '-') == 3)
          {
            mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 2)
            if (str_count(mirnaname, '-') == 3)
            {
              mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 1)
              if (str_count(mirnaname, '-') == 3)
              {
                mirnaname <- str_sub(mirnaname, 0, str_length(mirnaname) - 1)
                
              }
            }
            
          }
          mrna.search.results <-
            subset [grep((mirnaname), subset$probName, fixed = TRUE), ]
          
          length(mrna.search.results[, 1])
          # bind searched results
          if (nrow(mrna.search.results) > 0)
          {
            for (k in 1:nrow(mrna.search.results))
            {
              if ((mrna.search.results$probName[k]) == mirnaname)
              {
                sampleRow <- sample
                sampleRow$ProbName <- mrna.search.results$probName[k]
                sampleRow$Ftest <- mrna.search.results$Ftest[k]
                sampleRow$Ttest <- mrna.search.results$Ttest[k]
                sampleRow$Diff <- mrna.search.results$Diff[k]
                sampleRow$P_Val <- mrna.search.results$P_Val[k]
                sampleRow$P_Val_Less <-
                  mrna.search.results$P_Val_Less[k]
                sampleRow$P_Val_Greater <-
                  mrna.search.results$P_Val_Greater[k]
                
                result.part <- rbind(result.part, sampleRow)
                
              }
            }
            
            
          }
          
          
        }
      }
      
      
    }
    
    result.part
    
  }

# write.csv(file="~/Documents/Melonoma_TCGA-SKCM/2019_June/ttest.results_mirna_mrna.csv",r)

# result <- read.csv(file="~/Documents/Melonoma_TCGA-SKCM/2019_June/ttest.results_mirna_mrna.csv")

flagCrossTableOfSignificancyResults <- function(result) {
  result$prob_p[1] <- NA
  result$target_p[1] <- NA
  
  result$relation_type <-
    factor(result$relation_type , levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
  
  for (i in 1:nrow(result))
  {
    cat("\014")
    cat(i)
    if (!is.na(result$ProbName[i]))
    {
      if (result$Diff[i] == 1)
      {
        if (result$P_Val_Less[i] > result$P_Val_Greater[i])
        {
          result$prob_p[i] <- "L"
        } else
        {
          result$prob_p[i] <- "H"
        }
      }
      
    }
    
    if (!is.na(result$targetName[i]))
    {
      if (result$targetDiff[i] == 1)
      {
        if (result$targetP_Val_Less[i] > result$targetP_Val_Greater[i])
        {
          result$target_p[i] <- "L"
        } else
        {
          result$target_p[i] <- "H"
        }
      }
    }
    
    
    # TYPE-1 : miRNA H mRNA H
    
    if (result$prob_p[i] %in% "H"  &
        result$target_p[i] %in%   "H") {
      result$relation_type[i] <- 1
    }
    
    # TYPE-2 : miRNA H mRNA L
    if (result$prob_p[i] %in%  "H"  &
        result$target_p[i] %in% "L") {
      result$relation_type[i] <- 2
    }
    
    
    # TYPE-3 : miRNA NA mRNA L
    if (is.na(result$prob_p[i])  &
        result$target_p[i] %in%  "L") {
      result$relation_type[i] <- 3
    }
    
    
    # TYPE-4 : miRNA HNA mRNA H
    if (is.na(result$prob_p[i])  &
        result$target_p[i] %in%  "H") {
      result$relation_type[i] <- 4
    }
    
    
    # TYPE-5 : miRNA L mRNA H
    
    if (result$prob_p[i] %in%  "L"  &
        result$target_p[i] %in%  "H") {
      result$relation_type[i] <- 5
    }
    
    # TYPE-6 : miRNA L mRNA L
    
    if (result$prob_p[i] %in%  "L"  &
        result$target_p[i] %in%  "L") {
      result$relation_type[i] <- 6
    }
    
    # TYPE-7 : miRNA L mRNA NA
    
    if (result$prob_p[i] %in%  "L"  &
        is.na(result$target_p[i])) {
      result$relation_type[i] <- 7
    }
    
    
    # TYPE-8 : miRNA H mRNA NA
    if (result$prob_p[i] %in%  "H"  &&
        is.na(result$target_p[i])) {
      result$relation_type[i] <- 8
    }
    
    
    # TYPE-9 : miRNA NA mRNA NA
    if (is.na(result$prob_p[i])  &&
        is.na(result$target_p[i])) {
      result$relation_type[i] <- 9
    }
    
  }
  
  result
  
  
}



# for(i in 1608:nrow(significancyTestResults)) #nrol(significancyTestResults))
# {
#   cat("\014")
#   cat(paste("..", i, ".."))
#
#   # create mrna Part
#   phrase<-significancyTestResults$probName[i]
#
#   phrase <- str_sub(phrase, 0,str_locate(phrase,"\\.")[1]-1)
#
#   sample$targetName<-significancyTestResults$probName[i]
#   sample$targetFtest<-significancyTestResults$Ftest[i]
#   sample$targetTtest<-significancyTestResults$Ttest[i]
#   sample$targetDiff<-significancyTestResults$Diff[i]
#   sample$targetP_Val<-significancyTestResults$P_Val[i]
#   sample$targetP_Val_Less<-significancyTestResults$P_Val_Less[i]
#   sample$targetP_Val_Greater<-significancyTestResults$P_Val_Greater[i]
#   result<-rbind(result,sample)
#   #search target mRNA of given mirna
#   search.results<- na.omit( mrna_mirna_mapping [grep( phrase,mrna_mirna_mapping$mRNA ),])
#
#   if ( nrow(search.results)  >0)
#   {
#     for( y in 1:nrow(search.results) )
#     {
#
#       cat("\014")
#       cat(paste("..",i,"..", nrow(search.results) ,"/",y))
#       # for each  mrna search if there are any significacy result
#       mirnaname<-str_to_lower(search.results$miRNA[y])
#       if ( str_count(mirnaname,'-') ==3)
#       {
#         mirnaname <- str_sub(mirnaname, 0,str_length(mirnaname)-2)
#         if ( str_count(mirnaname,'-') ==3)
#         {
#           mirnaname <- str_sub(mirnaname, 0,str_length(mirnaname)-1)
#           if ( str_count(mirnaname,'-') ==3)
#           {
#             mirnaname <- str_sub(mirnaname, 0,str_length(mirnaname)-1)
#
#           }
#         }
#
#       }
#       mrna.search.results<- subset [grep((mirnaname),subset$probName, fixed = TRUE ),]
#
#       length(mrna.search.results[,1])
#       # bind searched results
#       if (nrow(mrna.search.results)>0)
#       {
#         for(k in 1:nrow(mrna.search.results))
#         {
#           if ( (mrna.search.results$probName[k]) == mirnaname)
#           {
#             sampleRow<-sample
#             sampleRow$ProbName <- mrna.search.results$probName[k]
#             sampleRow$Ftest <- mrna.search.results$Ftest[k]
#             sampleRow$Ttest <- mrna.search.results$Ttest[k]
#             sampleRow$Diff <- mrna.search.results$Diff[k]
#             sampleRow$P_Val <- mrna.search.results$P_Val[k]
#             sampleRow$P_Val_Less <- mrna.search.results$P_Val_Less[k]
#             sampleRow$P_Val_Greater <- mrna.search.results$P_Val_Greater[k]
#             result<-rbind(result,sampleRow)
#
#           }
#         }
#
#
#       }
#
#     }
#   }
#
#
# }

# subset<-significancyTestResults[significancyTestResults$Diff==1, ]


logToFile("Main", "UndersampleSecondlevelresults_all.csv", Sonuc_Matrix_undersample$results, "L" )
save(Sonuc_Matrix_undersample, file="Sonuc_Matrix_undersample_2ndLevel.RData")
?save

Sonuc_Matrix_undersample_p1 <- Sonuc_Matrix_undersample

Type_1_2 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_2", ]  
Type_1_3 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_3", ]  
Type_1_4 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_4", ]  
Type_1_5 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_5", ]  
Type_1_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_6", ]  
Type_1_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_7", ]  
Type_1_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1_8", ]  



output <- rbind(output, as.data.frame(  cbind("Type 1-2", CI(Type_1_2$Sensitivity, ci=0.95)[2], CI(Type_1_2$Specificity, ci=0.95)[2], CI(Type_1_2$AccuracyPValue, ci=0.95)[2],CI(Type_1_2$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 1-3", CI(Type_1_3$Sensitivity, ci=0.95)[2], CI(Type_1_3$Specificity, ci=0.95)[2], CI(Type_1_3$AccuracyPValue, ci=0.95)[2],CI(Type_1_3$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 1-4", CI(Type_1_4$Sensitivity, ci=0.95)[2], CI(Type_1_4$Specificity, ci=0.95)[2], CI(Type_1_4$AccuracyPValue, ci=0.95)[2],CI(Type_1_4$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 1-5", CI(Type_1_5$Sensitivity, ci=0.95)[2], CI(Type_1_5$Specificity, ci=0.95)[2], CI(Type_1_5$AccuracyPValue, ci=0.95)[2],CI(Type_1_5$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 1-6", CI(Type_1_6$Sensitivity, ci=0.95)[2], CI(Type_1_6$Specificity, ci=0.95)[2], CI(Type_1_6$AccuracyPValue, ci=0.95)[2],CI(Type_1_6$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 1-7", CI(Type_1_7$Sensitivity, ci=0.95)[2], CI(Type_1_7$Specificity, ci=0.95)[2], CI(Type_1_7$AccuracyPValue, ci=0.95)[2],CI(Type_1_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 1-8", CI(Type_1_8$Sensitivity, ci=0.95)[2], CI(Type_1_8$Specificity, ci=0.95)[2], CI(Type_1_8$AccuracyPValue, ci=0.95)[2],CI(Type_1_8$Accuracy, ci=0.95)[2])))




Type_2_3 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_3", ]  
Type_2_4 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_4", ]  
Type_2_5 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_5", ]  
Type_2_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_6", ]  
Type_2_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_7", ]  
Type_2_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_8", ]  

output <- rbind(output, as.data.frame(  cbind("Type 2-3", CI(Type_2_3$Sensitivity, ci=0.95)[2], CI(Type_2_3$Specificity, ci=0.95)[2], CI(Type_2_3$AccuracyPValue, ci=0.95)[2],CI(Type_2_3$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 2-4", CI(Type_2_4$Sensitivity, ci=0.95)[2], CI(Type_2_4$Specificity, ci=0.95)[2], CI(Type_2_4$AccuracyPValue, ci=0.95)[2],CI(Type_2_4$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 2-5", CI(Type_2_5$Sensitivity, ci=0.95)[2], CI(Type_2_5$Specificity, ci=0.95)[2], CI(Type_2_5$AccuracyPValue, ci=0.95)[2],CI(Type_2_5$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 2-6", CI(Type_2_6$Sensitivity, ci=0.95)[2], CI(Type_2_6$Specificity, ci=0.95)[2], CI(Type_2_6$AccuracyPValue, ci=0.95)[2],CI(Type_2_6$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 2-7", CI(Type_2_7$Sensitivity, ci=0.95)[2], CI(Type_2_7$Specificity, ci=0.95)[2], CI(Type_2_7$AccuracyPValue, ci=0.95)[2],CI(Type_2_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 2-8", CI(Type_2_8$Sensitivity, ci=0.95)[2], CI(Type_2_8$Specificity, ci=0.95)[2], CI(Type_2_8$AccuracyPValue, ci=0.95)[2],CI(Type_2_8$Accuracy, ci=0.95)[2])))


Type_3_4 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_4", ]  
Type_3_5 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_5", ]  
Type_3_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_6", ]  
Type_3_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_7", ]  
Type_3_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_8", ]  

output <- rbind(output, as.data.frame(  cbind("Type 3-4", CI(Type_3_4$Sensitivity, ci=0.95)[2], CI(Type_3_4$Specificity, ci=0.95)[2], CI(Type_3_4$AccuracyPValue, ci=0.95)[2],CI(Type_3_4$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 3-5", CI(Type_3_5$Sensitivity, ci=0.95)[2], CI(Type_3_5$Specificity, ci=0.95)[2], CI(Type_3_5$AccuracyPValue, ci=0.95)[2],CI(Type_3_5$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 3-6", CI(Type_3_6$Sensitivity, ci=0.95)[2], CI(Type_3_6$Specificity, ci=0.95)[2], CI(Type_3_6$AccuracyPValue, ci=0.95)[2],CI(Type_3_6$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 3-7", CI(Type_3_7$Sensitivity, ci=0.95)[2], CI(Type_3_7$Specificity, ci=0.95)[2], CI(Type_3_7$AccuracyPValue, ci=0.95)[2],CI(Type_3_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 3-8", CI(Type_3_8$Sensitivity, ci=0.95)[2], CI(Type_3_8$Specificity, ci=0.95)[2], CI(Type_3_8$AccuracyPValue, ci=0.95)[2],CI(Type_3_8$Accuracy, ci=0.95)[2])))

Type_4_5 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_4_5", ]  
Type_4_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_4_6", ]  
Type_4_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_4_7", ]  
Type_4_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_4_8", ]  

output <- rbind(output, as.data.frame(  cbind("Type 4-5", CI(Type_4_5$Sensitivity, ci=0.95)[2], CI(Type_4_5$Specificity, ci=0.95)[2], CI(Type_4_5$AccuracyPValue, ci=0.95)[2],CI(Type_4_5$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 4-6", CI(Type_4_6$Sensitivity, ci=0.95)[2], CI(Type_4_6$Specificity, ci=0.95)[2], CI(Type_4_6$AccuracyPValue, ci=0.95)[2],CI(Type_4_6$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 4-7", CI(Type_4_7$Sensitivity, ci=0.95)[2], CI(Type_4_7$Specificity, ci=0.95)[2], CI(Type_4_7$AccuracyPValue, ci=0.95)[2],CI(Type_4_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 4-8", CI(Type_4_8$Sensitivity, ci=0.95)[2], CI(Type_4_8$Specificity, ci=0.95)[2], CI(Type_4_8$AccuracyPValue, ci=0.95)[2],CI(Type_4_8$Accuracy, ci=0.95)[2])))

Type_5_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_5_6", ]  
Type_5_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_5_7", ]  
Type_5_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_5_8", ]  

output <- rbind(output, as.data.frame(  cbind("Type 5-6", CI(Type_5_6$Sensitivity, ci=0.95)[2], CI(Type_5_6$Specificity, ci=0.95)[2], CI(Type_5_6$AccuracyPValue, ci=0.95)[2],CI(Type_5_6$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 5-7", CI(Type_5_7$Sensitivity, ci=0.95)[2], CI(Type_5_7$Specificity, ci=0.95)[2], CI(Type_5_7$AccuracyPValue, ci=0.95)[2],CI(Type_5_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 5-8", CI(Type_5_8$Sensitivity, ci=0.95)[2], CI(Type_5_8$Specificity, ci=0.95)[2], CI(Type_5_8$AccuracyPValue, ci=0.95)[2],CI(Type_5_8$Accuracy, ci=0.95)[2])))


Type_6_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_6_7", ]  
Type_6_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_6_8", ]  

output <- rbind(output, as.data.frame(  cbind("Type 6-7", CI(Type_6_7$Sensitivity, ci=0.95)[2], CI(Type_6_7$Specificity, ci=0.95)[2], CI(Type_6_7$AccuracyPValue, ci=0.95)[2],CI(Type_6_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 6-8", CI(Type_6_8$Sensitivity, ci=0.95)[2], CI(Type_6_8$Specificity, ci=0.95)[2], CI(Type_6_8$AccuracyPValue, ci=0.95)[2],CI(Type_6_8$Accuracy, ci=0.95)[2])))


Type_7_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_7_8", ]  

output <- rbind(output, as.data.frame(  cbind("Type 7-8", CI(Type_6_8$Sensitivity, ci=0.95)[2], CI(Type_7_8$Specificity, ci=0.95)[2], CI(Type_7_8$AccuracyPValue, ci=0.95)[2],CI(Type_7_8$Accuracy, ci=0.95)[2])))

View(output)
# colnames(output) <- c("Definition","Sensitivity", "Specificity", "PValue","Accuracy" )

logToFile("Analysis", "UndersampleResults.csv", output, "L" )

par(mfrow = c(1, 2)) 
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_2$Sensitivity, Type_2$Specificity,
        Type_3$Sensitivity, Type_3$Specificity ,Type_4$Sensitivity, Type_4$Specificity,
        Type_5$Sensitivity, Type_5$Specificity ,Type_6$Sensitivity, Type_6$Specificity,
        Type_7$Sensitivity, Type_7$Specificity ,Type_8$Sensitivity, Type_8$Specificity,
        main=" Model Comparision" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c(" ", "Type 1"," ", "Type 2",
                " ", "Type 3"," ", "Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))    
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_1_2$Sensitivity, Type_1_2$Specificity ,
        Type_1_3$Sensitivity, Type_1_3$Specificity ,Type_1_4$Sensitivity, Type_1_4$Specificity,
        Type_1_5$Sensitivity, Type_1_5$Specificity ,Type_1_6$Sensitivity, Type_1_6$Specificity,
        Type_1_7$Sensitivity, Type_1_7$Specificity ,Type_1_8$Sensitivity, Type_1_8$Specificity,
        main=" Model Comparision for Type 1" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c( " ", "Type 1", " ", "Type 2",
                 " ", "Type 3"," ", "Type 4",
                 " ", "Type 5"," ", "Type 6 ",
                 " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))    



par(mfrow = c(1, 2)) 
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_2$Sensitivity, Type_2$Specificity,
        Type_3$Sensitivity, Type_3$Specificity ,Type_4$Sensitivity, Type_4$Specificity,
        Type_5$Sensitivity, Type_5$Specificity ,Type_6$Sensitivity, Type_6$Specificity,
        Type_7$Sensitivity, Type_7$Specificity ,Type_8$Sensitivity, Type_8$Specificity,
        main=" Model Comparision" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c(" ", "Type 1"," ", "Type 2",
                " ", "Type 3"," ", "Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))  
boxplot(Type_2$Sensitivity, Type_2$Specificity ,Type_2_3$Sensitivity, Type_2_3$Specificity ,Type_2_4$Sensitivity, Type_2_4$Specificity,
        Type_2_5$Sensitivity, Type_2_5$Specificity ,Type_2_6$Sensitivity, Type_2_6$Specificity,
        Type_2_7$Sensitivity, Type_2_7$Specificity ,Type_2_8$Sensitivity, Type_2_8$Specificity,
        main=" Model Comparision for Type 2" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17, 19,20),
        names=c(" ", "Type 2",
                 " ", "Type 3"," ", "Type 4",
                 " ", "Type 5"," ", "Type 6 ",
                 " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))    



par(mfrow = c(1, 2)) 
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_2$Sensitivity, Type_2$Specificity,
        Type_3$Sensitivity, Type_3$Specificity ,Type_4$Sensitivity, Type_4$Specificity,
        Type_5$Sensitivity, Type_5$Specificity ,Type_6$Sensitivity, Type_6$Specificity,
        Type_7$Sensitivity, Type_7$Specificity ,Type_8$Sensitivity, Type_8$Specificity,
        main=" Model Comparision" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c(" ", "Type 1"," ", "Type 2",
                " ", "Type 3"," ", "Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))  
boxplot(Type_3$Sensitivity, Type_3$Specificity, Type_3_4$Sensitivity, Type_3_4$Specificity,
        Type_3_5$Sensitivity, Type_3_5$Specificity ,Type_3_6$Sensitivity, Type_3_6$Specificity,
        Type_3_7$Sensitivity, Type_3_7$Specificity ,Type_3_8$Sensitivity, Type_3_8$Specificity,
        main=" Model Comparision for Type 3" ,
        at= c(1,2,4,5,7,8,10,11,13,14, 16,17),
        names=c(" ","Type 3", " ","Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))    


par(mfrow = c(1, 2)) 
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_2$Sensitivity, Type_2$Specificity,
        Type_3$Sensitivity, Type_3$Specificity ,Type_4$Sensitivity, Type_4$Specificity,
        Type_5$Sensitivity, Type_5$Specificity ,Type_6$Sensitivity, Type_6$Specificity,
        Type_7$Sensitivity, Type_7$Specificity ,Type_8$Sensitivity, Type_8$Specificity,
        main=" Model Comparision" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c(" ", "Type 1"," ", "Type 2",
                " ", "Type 3"," ", "Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))  
boxplot(Type_4$Sensitivity, Type_4$Specificity,Type_4_5$Sensitivity, Type_4_5$Specificity ,Type_4_6$Sensitivity, Type_4_6$Specificity,
        Type_4_7$Sensitivity, Type_4_7$Specificity ,Type_4_8$Sensitivity, Type_4_8$Specificity,
        main=" Model Comparision for Type 4" ,
        at= c(1,2,4,5,7,8,10,11,13,14),
        names=c( " ", "Type 4", " ", "Type 5"," ", "Type 6 ",
                 " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))    



par(mfrow = c(1, 2)) 
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_2$Sensitivity, Type_2$Specificity,
        Type_3$Sensitivity, Type_3$Specificity ,Type_4$Sensitivity, Type_4$Specificity,
        Type_5$Sensitivity, Type_5$Specificity ,Type_6$Sensitivity, Type_6$Specificity,
        Type_7$Sensitivity, Type_7$Specificity ,Type_8$Sensitivity, Type_8$Specificity,
        main=" Model Comparision" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c(" ", "Type 1"," ", "Type 2",
                " ", "Type 3"," ", "Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))  
boxplot(Type_5$Sensitivity, Type_5$Specificity,
        Type_5_6$Sensitivity, Type_5_6$Specificity,
        Type_5_7$Sensitivity,Type_5_7$Specificity ,
        Type_5_8$Sensitivity, Type_5_8$Specificity,
        main=" Model Comparision for Type 5" ,
        at= c(1,2,4,5,7,8,10,11),
        names=c( " ", "Type 5 ", " ", "Type 6 ",
                 " ", "Type 7" , " ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))   


par(mfrow = c(1, 2)) 
boxplot(Type_1$Sensitivity, Type_1$Specificity ,Type_2$Sensitivity, Type_2$Specificity,
        Type_3$Sensitivity, Type_3$Specificity ,Type_4$Sensitivity, Type_4$Specificity,
        Type_5$Sensitivity, Type_5$Specificity ,Type_6$Sensitivity, Type_6$Specificity,
        Type_7$Sensitivity, Type_7$Specificity ,Type_8$Sensitivity, Type_8$Specificity,
        main=" Model Comparision" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),
        names=c(" ", "Type 1"," ", "Type 2",
                " ", "Type 3"," ", "Type 4",
                " ", "Type 5"," ", "Type 6 ",
                " ", "Type 7"," ", "Type 8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))  
boxplot(Type_6$Sensitivity, Type_6$Specificity,Type_6_7$Sensitivity, Type_6_7$Specificity ,
        Type_6_8$Sensitivity, Type_6_8$Specificity,Type_7$Sensitivity, Type_7$Specificity,
        Type_8$Sensitivity, Type_8$Specificity,Type_7_8$Sensitivity, Type_7_8$Specificity ,
        main=" Model Comparision for Type 6/7/8" ,
        at= c(1,2,4,5,7,8,10,11,13,14,16,17),
        names=c(" ", "Type 6"," ", "Type 6-7",
                " ", "Type 6-8"," ", "Type 7",
                " ", "Type 8", " ", "Type 7-8"),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))    



### T TEST for all pairs...



t.test(Type_2$Accuracy, Type_2_6$Accuracy,alternative = c("two.sided", "less", "greater"), paired = FALSE, conf.level = 0.95)
t.test(Type_6$Accuracy, Type_2_6$Accuracy,alternative = c("two.sided", "less", "greater"), paired = FALSE, conf.level = 0.95)

result.list <- c(Type_1_2,Type_1_3)

vals <- NULL

getStats <- function (variable)
{
        Delta_M  <- abs(CI(variable$Sensitivity, ci=0.95)[2]- CI(variable$Specificity, ci=0.95)[2])
        Delta_A <- abs(CI(variable$Accuracy, ci=0.95)[1]-CI(variable$Accuracy, ci=0.95)[3])
        Delta_Sen <- abs(min(variable$Sensitivity)-max(variable$Sensitivity))
        Delta_Spe <- abs(min(variable$Specificity)-max(variable$Specificity))
        
        Acc <- CI(variable$Accuracy, ci=0.95)[2]
        Sen <- CI(variable$Sensitivity, ci=0.95)[2]
        Spe <- CI(variable$Specificity, ci=0.95)[2]
        vals <- cbind(variable$definition.type[1],Acc,Sen, Spe,Delta_A, Delta_M ,Delta_Sen,Delta_Spe  )      
        vals
}


vals <- rbind(vals,getStats(Type_1))
vals <- rbind(vals,getStats(Type_2))
vals <- rbind(vals,getStats(Type_3))    
vals <- rbind(vals,getStats(Type_4))
vals <- rbind(vals,getStats(Type_5))
vals <- rbind(vals,getStats(Type_6))    
vals <- rbind(vals,getStats(Type_7))
vals <- rbind(vals,getStats(Type_8))
    
vals <- rbind(vals,getStats(Type_1_2))
vals <- rbind(vals,getStats(Type_1_3))    
vals <- rbind(vals,getStats(Type_1_4))
vals <- rbind(vals,getStats(Type_1_5))
vals <- rbind(vals,getStats(Type_1_6))    
vals <- rbind(vals,getStats(Type_1_7))
vals <- rbind(vals,getStats(Type_1_8))

vals <- rbind(vals,getStats(Type_2_3))    
vals <- rbind(vals,getStats(Type_2_4))
vals <- rbind(vals,getStats(Type_2_5))
vals <- rbind(vals,getStats(Type_2_6))    
vals <- rbind(vals,getStats(Type_2_7))
vals <- rbind(vals,getStats(Type_2_8))

vals <- rbind(vals,getStats(Type_3_4))
vals <- rbind(vals,getStats(Type_3_5))
vals <- rbind(vals,getStats(Type_3_6))    
vals <- rbind(vals,getStats(Type_3_7))
vals <- rbind(vals,getStats(Type_3_8))

vals <- rbind(vals,getStats(Type_4_5))
vals <- rbind(vals,getStats(Type_4_6))    
vals <- rbind(vals,getStats(Type_4_7))
vals <- rbind(vals,getStats(Type_4_8))

vals <- rbind(vals,getStats(Type_5_6))    
vals <- rbind(vals,getStats(Type_5_7))
vals <- rbind(vals,getStats(Type_5_8))

vals <- rbind(vals,getStats(Type_6_7))
vals <- rbind(vals,getStats(Type_6_8))

vals <- rbind(vals,getStats(Type_7_8))

logToFile("Analysis","ModelComparePreSets.csv", vals,"L")

View(vals)

f  <- Type_1_7
s  <- Type_3_7
n  <- "Type_1_7_3"

Type_XXX <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% n, ] 
Type_XXX <- Type_XXX[Type_XXX$`model_list[j]` %in% "svmLinear",]

CI(Type_XXX$AccuracyPValue, ci=0.95)[2]

t.test(f$Accuracy, Type_XXX$Accuracy,alternative = c("two.sided", "less", "greater"), paired = FALSE, conf.level = 0.95)
t.test(s$Accuracy, Type_XXX$Accuracy,alternative = c("two.sided", "less", "greater"), paired = FALSE, conf.level = 0.95)

getStats(Type_XXX)
getStats(f)
getStats(s)

par(mfrow = c(1, 1)) 
boxplot(Type_XXX$Sensitivity, Type_XXX$Specificity,
        f$Sensitivity, f$Specificity ,
        s$Sensitivity, s$Specificity,
        main=" Model Comparision for Type 1-7 & 3_7" ,
        at= c(1,2,4,5,7,8),
        names=c(" ", n,
                " ", " 3-7",
                " ", "6_7" ),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))   

getResults<- function (variable, Model)
{
        selectedSet <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% variable, ] 
        selectedSet <- selectedSet[selectedSet$`model_list[j]` %in% Model,]
        
        Delta_M  <- abs(CI(selectedSet$Sensitivity, ci=0.95)[2]- CI(selectedSet$Specificity, ci=0.95)[2])
        Delta_A <- abs(CI(selectedSet$Accuracy, ci=0.95)[1]-CI(selectedSet$Accuracy, ci=0.95)[3])
        Delta_Sen <- abs(min(selectedSet$Sensitivity)-max(selectedSet$Sensitivity))
        Delta_Spe <- abs(min(selectedSet$Specificity)-max(selectedSet$Specificity))
        
        Acc <- CI(selectedSet$Accuracy, ci=0.95)[2]
        Sen <- CI(selectedSet$Sensitivity, ci=0.95)[2]
        Spe <- CI(selectedSet$Specificity, ci=0.95)[2]
        pVals <- CI(selectedSet$AccuracyPValue, ci=0.95)[2]
        vals <- cbind(selectedSet$definition.type[1],Model, Acc,Sen, Spe,pVals, Delta_A, Delta_M ,Delta_Sen,Delta_Spe )      
        vals
}

mySubSets <- c("Type_1_7_3","Type_2_6_8","Type_2_6_7","Type_2_3_8","Type_1_7_8","Type_1_7_6","Type_1_7_2","Type_2_3_7")

results_ <- NULL

for (i in 1: length(mySubSets)) {
        results_ <- rbind(results,getResults(mySubSets[i],"svmLinear"))
}


mySubSets <- c("Type_3_7_8")

for (variable in mySubSets) {
        results_ <- rbind(results_,getResults(variable,"svmLinear"))
        results_ <- rbind(results_,getResults(variable,"svmRadial"))
        results_ <- rbind(results_,getResults(variable,"svmPoly"))
        results_ <- rbind(results_,getResults(variable,"nnet"))
        
}

View(results_ )

logToFile("Analysis", "3rdLevel.csv",results,"L")


############# Box Plot

Type_2_3_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_3_6", ] 
Type_2_3_6 <- Type_2_3_6[Type_2_3_6$`model_list[j]` %in% "svmLinear",]

Type_6_7_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_6_7_8", ] 
Type_6_7_8 <- Type_6_7_8[Type_6_7_8$`model_list[j]` %in% "svmLinear",]

Type_3_6_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_6_7", ] 
Type_3_6_7 <- Type_3_6_7[Type_3_6_7$`model_list[j]` %in% "svmLinear",]

Type_3_7_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3_7_8", ] 
Type_3_7_8 <- Type_3_7_8[Type_3_7_8$`model_list[j]` %in% "svmLinear",]

Type_2_3_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2_3_7", ] 
Type_2_3_7 <- Type_2_3_7[Type_2_3_7$`model_list[j]` %in% "svmLinear",]


boxplot(Type_2_3_6$Sensitivity, Type_2_3_6$Specificity,
        Type_6_7_8$Sensitivity, Type_6_7_8$Specificity ,
        Type_3_6_7$Sensitivity, Type_3_6_7$Specificity,
        Type_3_7_8$Sensitivity, Type_3_7_8$Specificity,
        Type_2_3_7$Sensitivity, Type_2_3_7$Specificity,
        main=" Model Comparision for selected pairs" ,
        at= c(1,2,4,5,7,8,10,11,13,14),
        names=c(" ", "2_3_6",
                " ", "6_7_8",
                " ", "3_6_7",
                " ", "3_7_8" ,
                " ", "2_3_7" ),
        col = c("orange","red"),
        las=2,
        border = "brown",
        horizontal = TRUE,
        notch = FALSE)
mtext(paste("Orange:Sensitivity   Red: Specificity"))
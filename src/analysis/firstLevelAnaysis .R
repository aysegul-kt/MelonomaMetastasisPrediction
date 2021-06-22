


logToFile("Main", "UndersampleFirstlevelresults_all.csv", Sonuc_Matrix_undersample$results, "L" )
save(Sonuc_Matrix_undersample, file="Sonuc_Matrix_undersample_1stLevel.RData")
?save

Sonuc_Matrix_undersample_1stlevel <- Sonuc_Matrix_undersample

Type_1 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_1", ]  
Type_2 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_2", ]  
Type_3 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_3", ]  
Type_4 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_4", ]  
Type_5 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_5", ]  
Type_6 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_6", ]  
Type_7 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_7", ]  
Type_8 <- Sonuc_Matrix_undersample$results[Sonuc_Matrix_undersample$results$definition.type %in% "Type_8", ]  

op <- par(mfrow = c(1, 2)) 
par(op)
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


output <- as.data.frame(  cbind("Type 1", CI(Type_1$Sensitivity, ci=0.95)[2], CI(Type_1$Specificity, ci=0.95)[2], CI(Type_1$AccuracyPValue, ci=0.95)[2],CI(Type_1$Accuracy, ci=0.95)[2]))
output <- rbind(output, as.data.frame(  cbind("Type 2", CI(Type_2$Sensitivity, ci=0.95)[2], CI(Type_2$Specificity, ci=0.95)[2], CI(Type_2$AccuracyPValue, ci=0.95)[2],CI(Type_2$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 3", CI(Type_3$Sensitivity, ci=0.95)[2], CI(Type_3$Specificity, ci=0.95)[2], CI(Type_3$AccuracyPValue, ci=0.95)[2],CI(Type_3$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 4", CI(Type_4$Sensitivity, ci=0.95)[2], CI(Type_4$Specificity, ci=0.95)[2], CI(Type_4$AccuracyPValue, ci=0.95)[2],CI(Type_4$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 5", CI(Type_5$Sensitivity, ci=0.95)[2], CI(Type_5$Specificity, ci=0.95)[2], CI(Type_5$AccuracyPValue, ci=0.95)[2],CI(Type_5$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 6", CI(Type_6$Sensitivity, ci=0.95)[2], CI(Type_6$Specificity, ci=0.95)[2], CI(Type_6$AccuracyPValue, ci=0.95)[2],CI(Type_6$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 7", CI(Type_7$Sensitivity, ci=0.95)[2], CI(Type_7$Specificity, ci=0.95)[2], CI(Type_7$AccuracyPValue, ci=0.95)[2],CI(Type_7$Accuracy, ci=0.95)[2])))
output <- rbind(output, as.data.frame(  cbind("Type 8", CI(Type_8$Sensitivity, ci=0.95)[2], CI(Type_8$Specificity, ci=0.95)[2], CI(Type_8$AccuracyPValue, ci=0.95)[2],CI(Type_8$Accuracy, ci=0.95)[2])))



#colnames(output) <- c("Definition","Sensitivity", "Specificity", "PValue","Accuracy" )

logToFile("Main", "UndersampleFirstlevelresults.csv", output, "L" )

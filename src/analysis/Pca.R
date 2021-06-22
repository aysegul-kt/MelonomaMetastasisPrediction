library("FactoMineR")
library("factoextra")

#####################Descritptions#####################
# computes principal component analysis on the active individuals/variables:
res.pca <- PCA(decathlon2.active, graph = FALSE)
print(res.pca)


#  The eigenvalues and the proportion of variances (i.e., information) 
#retained by the principal components (PCs) can be extracted using the function
# get_eigenvalue() [factoextra package].
eig.val <- get_eigenvalue(res.pca)
eig.val

#The scree plot can be produced using the function fviz_eig() or fviz_screeplot() [factoextra package].
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#A simple method to extract the results, for variables, from a PCA output is to use the functionget_pca_var() 
#[factoextra package]. This function provides a list of matrices containing all
#the results for the active variables (coordinates, correlation between variables and axes, squared cosine and contributions)

var <- get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

# Coordinates of variables
head(var$coord, 4)

# To plot variables, type this:
fviz_pca_var(res.pca, col.var = "black")

# Quality of representation
# The quality of representation of the variables on factor map is called cos2 (square cosine, squared coordinates)

library("corrplot")
head(var$cos2, 4)
corrplot(var$cos2, is.corr=FALSE)


# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

#the function dimdesc() [in FactoMineR], for dimension description, can be used to 
#identify the most significantly associated variables with a given principal component . It can be used as follow:

res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1

res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1


#Color by groups

iris.pca <- PCA(iris[,-5], graph = FALSE)

fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


#############Application##################

lvariables <- rbind(pool$miRNAFeatures,pool$mRNAFeatures,pool$methylationFeatures)
setorder(lvariables,P_Val) 

lvariablesList <- as.vector(lvariables[1:2000,]$probName)  


getPCAImportance <- function(variables)
{

modelInpts <-
  fxGetModelDataInputs (
    dataSample = dataSampleList[[1]],
    mode = "vars",
    variableList = variables,
    runtypes = c("sm", "us"),
    usePrev = FALSE,
    data = NULL,
    preruuid = NULL
  )


#calculate PCA
data.pca <- PCA(modelInpts$sm$trainning[,-2001], graph = FALSE,ncp = 175)
fviz_eig(data.pca, addlabels = TRUE, ylim = c(0, 35), ncp = 30,choice = "variance",geom="bar")
head(data.pca$eig[1:20,])
plot(data.pca$eig[,2], xlab = "Principal Component", ylab = " Variance Explained", type="b")

plot(data.pca$eig[,3], xlab = "Principal Component (#)", ylab = "Cumulative Proportion of Variance Explained (%)", type="b")
abline(h=98,v=175) # 193 for step 2

pca.var <- get_pca_var(data.pca)

fviz_contrib(data.pca, choice = "var", axes = 120, top = 200)

head(pca.var$contrib[1:20,])
library("corrplot")

corrplot(pca.var$contrib, is.corr=FALSE, tl.cex = 0.5, method = "color") 


res.desc <- dimdesc(data.pca, axes = c(1:175), proba = 0.0001)

View(res.desc)
ContributionsFrame <- data.frame("attribute"=as.character(NULL),"dim"=as.double(NULL),"correlation"=as.double(NULL),"pvalue"=as.double(NULL))

for (i in 1:175) {
  if(!is.null(res.desc[[i]]$quanti)){
    ContributionsFrame <- rbind(ContributionsFrame,data.frame("attribute"=as.character( rownames(res.desc[[i]]$quanti)), "dim"=c(rep(i,nrow(res.desc[[i]]$quanti))), "correlation"=res.desc[[i]]$quanti[,1], "pvalue"=res.desc[[i]]$quanti[,2]))
    
  }
  
}

row.names(ContributionsFrame) <- NULL

View(ContributionsFrame)
ContributionsFrame

}





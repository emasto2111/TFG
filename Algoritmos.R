# Ingenieria de la salud, Emanuela Maria Stoia

# Library imports:
library(DESeq2)
library(rlog)
library(WGCNA)
library(MCL)
library(igraph)
library(dplyr)
library(factoextra)
library(cluster)
require(ppclust)
require(fclust)
library("dbscan")


# Directory:
setwd("~/TFG/CodigoTFG/CodigoRstudio/Data")

# Algortihms:

## WGCNA
toc_to_mod_WGCNA <- function(data,powers) {
  sft <- pickSoftThreshold(toc, powerVector = powers, verbose = 5)
  
  min_fit <- 0.9
  if(any(sft$fitIndices$SFT.R.sq > min_fit) == FALSE) { 
    warning(" no power values gave a 0.9 fit value therefore, 
            we are going to use the fit value of maximum power") 
    power <- 20
  } else {
    power <- min(which(sft$fitIndices$SFT.R.sq > min_fit))
  }
  
  
  
  data<- t(data)
  # WGCNA code to go from toc to gene-module assignations
  net = blockwiseModules(data, power = power,
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE, 
                         verbose = 3)
  datos_fin<- data.frame(colnames(data),net$colors)
  return(datos_fin)
  
}

## MCL
toc_to_mod_MCL <- function(data) {
  data<- scale(data)
  tdata <- t(data)
  matrixadj <- adjacency(tdata)
  modules <- mcl(x = matrixadj, addLoops = TRUE)
  nGenes <- rownames(data) 
  datos_fin <- data.frame(modules$Cluster, nGenes)
  return(datos_fin)
}

## FUZZY
toc_to_mod_FUZZY <- function(data, centers) {
  res.fcm <- fcm(data, centers = centers)
  nGenes <- rownames(data)
  datos_fin <- data.frame(nGenes,res.fcm$cluster)
  return(datos_fin)
}

## F4 K-MEDIOSIS

toc_to_mod_KMEDOIDS <- function(data) {
  kmed <- pam(data, k = 6)
  clust <- fviz_cluster(kmed, data = data)
  datos_fin <- data.frame(clust$data$name, clust$data$cluster)
  return(datos_fin)
}
## F5 DbScan

toc_to_mod_DbScan <- function(data,eps) {
  MinPts<- ncol(data) + 1
  cl<- dbscan(data,eps=eps,MinPts=MinPts)
  nGenes <- rownames(data)
  datos_fin <- data.frame(nGenes, cl$cluster)
  return(datos_fin)
}
## Calculate time of duration to clustering different partition of data 
time_systemFor_method <- function(data,i,func,pow){
  listt<- list()
  datalist <- list()
  for (j in 1:length(i)){
    d<- data[sample(nrow(data),i[j]),]
    datalist<- system.time(func(d,pow))
    listt[[j]]<- c(datalist[1],datalist[2],datalist[3])
  }
  return(listt)
}
time_systemFor_methodX <- function(data,i,func){
  listt<- list()
  datalist <- list()
  for (j in 1:length(i)){
    d<- data[sample(nrow(data),i[j]),]
    datalist<- system.time(func(d))
    listt[[j]]<- c(datalist[1],datalist[2],datalist[3])
  }
  return(listt)
}

## Calculate memory to clustering different partition of data:

memory_systemFor_method <- function(data,i,func,pow){
  listt<- list()
  datalist <- list()
  for (j in 1:length(i)){
    d<- data[sample(nrow(data),i[j]),]
    a <- gc(reset = TRUE)
    datalist<- func(d,pow)
    b<- gc()
    listt[[j]]<- c(b[11]-a[11],b[12]-a[12])
  }
  return(listt)
}
memory_systemFor_methodX <- function(data,i,func){
  listt<- list()
  datalist <- list()
  for (j in 1:length(i)){
    d<- data[sample(nrow(data),i[j]),]
    a <- gc(reset = TRUE)
    datalist<- func(d)
    b<- gc()
    listt[[j]]<- c(b[11]-a[11],b[12]-a[12])
  }
  return(listt)
}


# Convert data into matrix 

input_file <- "counts.csv"
set.seed(1000)
toc <- read.csv(input_file, header=TRUE, sep= "")
toc <- toc[rowSums(toc)>20,]
toc2 <- toc[apply(toc, 1, sd)>10,]
condition <- factor(c(rep("control", each =4), rep("treated", each=3 )))
pe <- factor(c(rep("col0", each = 4),rep("col1", each = 3)))
targets <- data.frame(condition,pe)
dds=DESeqDataSetFromMatrix(toc2, targets, design=~condition)
toc <- rlog(dds)
toc<- assay(toc)


# Execute code with sample of 2000 gene
tocS<- toc[sample(nrow(toc), 2000),]
pow = c(c(1:10), seq(from = 12, to=20, by=1))
geneWGCNA <- toc_to_mod_WGCNA(tocS,pow )
geneMCL <- toc_to_mod_MCL(tocS)
geneFuzzy <- toc_to_mod_FUZZY(tocS,3)
geneKMed <- toc_to_mod_KMEDOIDS(tocS)
eps<- 0.85
geneDbScan <- toc_to_mod_DbScan(tocS,eps)

# Save data
setwd("~/TFG/CodigoTFG/Recuentos")
write.csv(geneWGCNA, file="geneWGCNA2000")
write.csv(geneMCL, file="geneMCL2000")
write.csv(geneFuzzy, file="geneFuzzy2000")
write.csv(geneKMed, file="geneKMedoids2000")
write.csv(geneDbScan, file="geneDbScan2000")

# Execute code with all of genes
pow = c(c(1:10), seq(from = 12, to=20, by=1))
geneWGCNA <- toc_to_mod_WGCNA(toc,pow )
#geneMCL <- toc_to_mod_MCL(toc)
geneFuzzy <- toc_to_mod_FUZZY(toc,4)
geneKMedoids <- toc_to_mod_KMEDOIDS(toc)
eps<- 0.35
geneDbScan <- toc_to_mod_DbScan(toc,eps)

# Save data
setwd("~/TFG/CodigoTFG/Recuentos")
write.csv(geneWGCNA, file="geneWGCNAall")
write.csv(geneFuzzy, file="geneFUZZYall")
write.csv(geneKMedoids, file="geneKMedoidsall")
write.csv(geneDbScan, file="geneDbScanall")

# Execution time
i<-c(1000,5000,7000,dim(toc)[1])
j<- c(1000,5000,7000)
timeMCL<- time_systemFor_methodX(toc,j,toc_to_mod_MCL)
timeWGCNA <- time_systemFor_method(toc,i,toc_to_mod_WGCNA,pow)
timeFUZZY <- time_systemFor_method(toc,i,toc_to_mod_FUZZY,3)
timeKMEDOIDS <- time_systemFor_methodX(toc,i,toc_to_mod_KMEDOIDS)
timeDBSCAN <- time_systemFor_method(toc,i,toc_to_mod_DbScan,0.65)

Split1000User<- c(timeWGCNA[[1]][1],timeFUZZY[[1]][1],timeKMEDOIDS[[1]][1],timeDBSCAN[[1]][1], timeMCL[[1]][1])
Split1000Elapsed<- c(timeWGCNA[[1]][3],timeFUZZY[[1]][3],timeKMEDOIDS[[1]][3],timeDBSCAN[[1]][3], timeMCL[[1]][3])
Split5000User<- c(timeWGCNA[[2]][1],timeFUZZY[[2]][1],timeKMEDOIDS[[2]][1],timeDBSCAN[[2]][1], timeMCL[[2]][1])
Split5000Elapsed<- c(timeWGCNA[[2]][3],timeFUZZY[[2]][3],timeKMEDOIDS[[2]][3],timeDBSCAN[[2]][3], timeMCL[[2]][3])
Split7000User<- c(timeWGCNA[[3]][1],timeFUZZY[[3]][1],timeKMEDOIDS[[3]][1],timeDBSCAN[[3]][1], timeMCL[[3]][1])
Split7000Elapsed<- c(timeWGCNA[[3]][3],timeFUZZY[[3]][3],timeKMEDOIDS[[3]][3],timeDBSCAN[[3]][3], timeMCL[[3]][3])
SplitDimUser<- c(timeWGCNA[[4]][1],timeFUZZY[[4]][1],timeKMEDOIDS[[4]][1],timeDBSCAN[[4]][1], NA)
SplitDimElapsed<- c(timeWGCNA[[4]][3],timeFUZZY[[4]][3],timeKMEDOIDS[[4]][3],timeDBSCAN[[4]][3], NA)
rw <- c("WGCNA", "FUZZY", "KMEDOIDS", "DBSCAN", "MCL")
cl <- c("1000User","1000Elapsed","5000User","5000Elapsed","7000User","7000Elapsed","DimUser","DimElapsed")
timedf <- data.frame(Split1000User,Split1000Elapsed,Split5000User,Split5000Elapsed,Split7000User,Split7000Elapsed,SplitDimUser,SplitDimElapsed)
rownames(timedf)<- rw
colnames(timedf)<- cl
timedf<- t(timedf)

setwd("~/TFG/CodigoTFG/Tiempos")
write.csv(timedf, file="timeExecution")

# Memory Usage
i<-c(1000,5000,7000,dim(toc)[1])
j<- c(1000,5000,7000)
memoryMCL<- memory_systemFor_methodX(toc,j,toc_to_mod_MCL)
memoryWGCNA <- memory_systemFor_method(toc,i,toc_to_mod_WGCNA,pow)
memoryFUZZY <- memory_systemFor_method(toc,i,toc_to_mod_FUZZY,3)
memoryKMEDOIDS <- memory_systemFor_methodX(toc,i,toc_to_mod_KMEDOIDS)
memoryDBSCAN <- memory_systemFor_method(toc,i,toc_to_mod_DbScan,0.65)

Split1000Ncells<- c(memoryWGCNA[[1]][1],memoryFUZZY[[1]][1],memoryKMEDOIDS[[1]][1],memoryDBSCAN[[1]][1], memoryMCL[[1]][1])
Split1000Vcells<- c(memoryWGCNA[[1]][2],memoryFUZZY[[1]][2],memoryKMEDOIDS[[1]][2],memoryDBSCAN[[1]][2], memoryMCL[[1]][2])
Split5000Ncells<- c(memoryWGCNA[[2]][1],memoryFUZZY[[2]][1],memoryKMEDOIDS[[2]][1],memoryDBSCAN[[2]][1], memoryMCL[[2]][1])
Split5000Vcells<- c(memoryWGCNA[[2]][2],memoryFUZZY[[2]][2],memoryKMEDOIDS[[2]][2],memoryDBSCAN[[2]][2], memoryMCL[[2]][2])
Split7000Ncells<- c(memoryWGCNA[[3]][1],memoryFUZZY[[3]][1],memoryKMEDOIDS[[3]][1],memoryDBSCAN[[3]][1], memoryMCL[[3]][1])
Split7000Vcells<- c(memoryWGCNA[[3]][2],memoryFUZZY[[3]][2],memoryKMEDOIDS[[3]][2],memoryDBSCAN[[3]][2], memoryMCL[[3]][2])
SplitDimNcells<- c(memoryWGCNA[[4]][1],memoryFUZZY[[4]][1],memoryKMEDOIDS[[4]][1],memoryDBSCAN[[4]][1], NA)
SplitDimVcells<- c(memoryWGCNA[[4]][2],memoryFUZZY[[4]][2],memoryKMEDOIDS[[4]][2],memoryDBSCAN[[4]][2], NA)

rw <- c("WGCNA", "FUZZY", "KMEDOIDS", "DBSCAN", "MCL")
cl <- c("1000NCells","1000VCells","5000NCells","5000VCells","7000NCells","7000VCells","DimNCells","DimVCells")
memoriadf <- data.frame(Split1000Ncells,Split1000Vcells,Split5000Ncells,Split5000Vcells,Split7000Ncells,Split7000Vcells,SplitDimNcells,SplitDimVcells)
rownames(memoriadf)<- rw
colnames(memoriadf)<- cl
memoriadf<- t(memoriadf)

setwd("~/TFG/CodigoTFG/Memoria")
write.csv(memoriadf, file="MemoryExecution")

## https://stackoverflow.com/questions/63361306/measure-peak-memory-usage-in-r

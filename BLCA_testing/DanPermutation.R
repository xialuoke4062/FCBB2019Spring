library(stringr)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

bestTable = read.csv("bigSummaryDf_BLCA1stTry.csv")
bestTable$X = NULL

newTable = bestTable[bestTable$co_occurance > 10, ]

OGtable = read.csv("updated_encoding_BLCA.csv")
rownames(OGtable) = OGtable$X
OGtable$X = NULL

permutateNull_sig <- function(tempMatrix) {
    numberOfCo = c()
    numberOfMu = c()
    
    for(i in 1:1000) {
        newTemp = tempMatrix
        newTemp[1, ] = sample(newTemp[1, ])
        newTemp = newTemp[, !(newTemp[1, ] == 0 & newTemp[2, ] == 0)] 
        
        c_occurance = sum(newTemp[1, ] == newTemp[2, ])
        m_exclusive = ncol(newTemp) - c_occurance
        
        numberOfCo = c(numberOfCo, c_occurance)
        numberOfMu = c(numberOfMu, m_exclusive)
    }
    
    returnArray = c(as.numeric(quantile(numberOfCo, probs = 0.95)), as.numeric(quantile(numberOfMu, probs = 0.95)))
    
    #return 
    returnArray
}

permutateNull_sigParallel <- function(tempMatrix) {
    numberOfCo = c()
    numberOfMu = c()
    parallelTest <- foreach(icount(1000), .combine=cbind) %dopar% {
        newTemp = tempMatrix
        newTemp[1, ] = sample(newTemp[1, ])
        newTemp = newTemp[, !(newTemp[1, ] == 0 & newTemp[2, ] == 0)] 
        
        c_occurance = sum(newTemp[1, ] == newTemp[2, ])
        c_occurance
    }
    
    parallelTest = as.vector(parallelTest[1, ])
    returnNum = as.numeric(quantile(parallelTest, probs = 0.95))
    #return 
    returnNum
}

minCo = c()
for(index in 1:nrow(newTable)){
    gene1 = as.character(newTable[index, "gene1"])
    gene2 = as.character(newTable[index, "gene2"])
    
    tempMatrix =OGtable[c(gene1, gene2), ]
    
    returnArray = permutateNull_sigParallel(tempMatrix)
    print(returnArray)
    minCo = c(minCo, returnArray)
    
}
newTable$cutOff = minCo

write.csv(newTable, file = "newTable_withCutoff.csv")

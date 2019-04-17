library(stringr)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#setwd("C:/Users/Dan/Dropbox (CahanLab)/Foundation_FinalProject")
setwd("~/Dropbox (CahanLab)/Foundation_FinalProject")
myTable = read.csv("updated_encoding_BLCA.csv")
rownames(myTable) = myTable$X
myTable$X = NULL

permutateNull_sigParallel <- function(tempMatrix) {
    numberOfCo = c()
    numberOfMu = c()
    parallelTest <- foreach(icount(500), .combine=cbind) %dopar% {
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


binomial_conversion <- function(myTable) {
    driverGenes = rownames(myTable)
    
    bigSummaryDf = data.frame(gene1 = NULL, gene2 = NULL, co_occurance = NULL, mutual_exclusion = NULL, co_cut = NULL)
    
    for(index1 in 1:length(driverGenes)) {
        for(index2 in 1:length(driverGenes)) {
            if (index2 > index1 ) {
                gene1 = driverGenes[index1]
                gene2 = driverGenes[index2]
                
                tempMatrix = myTable[c(gene1, gene2), ]
                
                # filter out incidence where everything is 0 
                tempMatrix = tempMatrix[, !(tempMatrix[1, ] == 0 & tempMatrix[2, ] == 0)] 
                
                c_occurance = sum(tempMatrix[1, ] == tempMatrix[2, ])
                m_exclusive = ncol(tempMatrix) - c_occurance

                sigCutoff = permutateNull_sigParallel(tempMatrix)
                
                fillerDf = data.frame(gene1 = gene1, gene2 = gene2, co_occurance = c_occurance, mutual_exclusion = m_exclusive, co_cut = sigCutoff)
                bigSummaryDf = rbind(bigSummaryDf, fillerDf)
                
                print(fillerDf)
            }
            
        }
    }

    #return 
    bigSummaryDf
}

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

mySummary = binomial_conversion(myTable)
write.csv(mySummary, file = "bigSummaryDf_BLCA.csv")


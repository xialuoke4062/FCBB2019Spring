# the samples that I need are 
library(stringr)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


# SKCM - Done 
# SARC - DONE 
# PAAD - DONE
# LUSC
# LUAD
# KIRP
# KIRC 
# CESC
# UCEC
# STAD 

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
   #returnNum = as.numeric(quantile(parallelTest, probs = 0.95))
   #return 
   parallelTest
   #returnNum
}

#cancerList = c("LGG", "GBM", "BRCA", "LIHC", "LAML", "SKCM", "SARC", "PAAD", "LUSC", "LUAD", "KIRP", "KIRC", "CESC","STAD", "UCEC", "BLCA") # first round finished up to LAML

#cancerList = c("KIRP", "KIRC", "CESC", "STAD", "SKCM", "SARC", "PAAD", "LUSC", "LUAD", "UCEC") 

cancerList = c("STAD", "SARC", "PAAD", "LUSC", "LUAD") 

for(myCancer in cancerList) {

    print(paste0("start ", myCancer))
    setwd(paste0("/Users/danpeng/Desktop/FCBB2019Spring/", myCancer, "/"))
    myTable = read.csv(paste0("encoded", myCancer, ".csv")) # change this to appropriate file name 

    rownames(myTable) = myTable$X
    myTable$X = NULL

    driverGenes = rownames(myTable)

    bigSummaryDf = data.frame(gene1 = NULL, gene2 = NULL, co_occurance = NULL, mutual_exclusion = NULL, co_cut = NULL, pval = NULL)

    distributionList = list()
    for(index1 in 1:length(driverGenes)) {
        for(index2 in 1:length(driverGenes)) {
            if (index2 > index1 ) {
                gene1 = driverGenes[index1]
                gene2 = driverGenes[index2]
                
                tempMatrix = myTable[c(gene1, gene2), ]
                # filter out incidence where everything is 0 
                if(sum(!(tempMatrix[1, ] == 0 & tempMatrix[2, ] == 0)) < 2) {
                    next
                }

                OGtempMatrix = tempMatrix

                tempMatrix = tempMatrix[, !(tempMatrix[1, ] == 0 & tempMatrix[2, ] == 0)] 
                
                c_occurance = sum(tempMatrix[1, ] == tempMatrix[2, ])
                m_exclusive = ncol(tempMatrix) - c_occurance
                
                if(myCancer == "STAD") {
                    minPairs = 2
                } else {
                    minPairs = 10
                }
                if(c_occurance <= minPairs) {
                    next
                }

                print(paste0("start ", myCancer))
                
                permuteSeq = permutateNull_sigParallel(OGtempMatrix)
                sigCutoff = as.numeric(quantile(permuteSeq, probs = 0.95))
                
                distributionList[[paste0(gene1, "_", gene2)]] = permuteSeq
                
                pval = sum(permuteSeq >= c_occurance) / length(permuteSeq)

                fillerDf = data.frame(gene1 = gene1, gene2 = gene2, co_occurance = c_occurance, mutual_exclusion = m_exclusive, co_cut = sigCutoff, pval = pval)
                
                bigSummaryDf = rbind(bigSummaryDf, fillerDf)

                write.csv(bigSummaryDf, file = paste0("bigSummaryDf_2ndtry", myCancer, ".csv"))
                save(distributionList, file = paste0("distributionList_2ndtry",myCancer,".rda"))
            }
        }
    }
}

library(stringr)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

setwd("/Users/danpeng/Desktop/FCBB2019Spring/BRCA/")
myTable = read.csv("encodedBRCA.csv") # change this to appropriate file name 
rownames(myTable) = myTable$X
myTable$X = NULL

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
   #returnNum = as.numeric(quantile(parallelTest, probs = 0.95))
   #return 
   parallelTest
   #returnNum
}

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
         
         tempMatrix = tempMatrix[, !(tempMatrix[1, ] == 0 & tempMatrix[2, ] == 0)] 
         
         c_occurance = sum(tempMatrix[1, ] == tempMatrix[2, ])
         m_exclusive = ncol(tempMatrix) - c_occurance
         
         if(c_occurance < 10) {
            next
         }
         
         #sigCutoff = permutateNull_sigParallel(tempMatrix)
         permuteSeq = permutateNull_sigParallel(tempMatrix)
         sigCutoff = as.numeric(quantile(permuteSeq, probs = 0.95))
         
         distributionList[[paste0(gene1, "_", gene2)]] = permuteSeq
         
         pval = sum(permuteSeq >= c_occurance) / length(permuteSeq)
         
         print(paste0(sum(permuteSeq >= c_occurance), " numbers of corruances"))
         print(length(permuteSeq))
         fillerDf = data.frame(gene1 = gene1, gene2 = gene2, co_occurance = c_occurance, mutual_exclusion = m_exclusive, co_cut = sigCutoff, pval = pval)
         print(fillerDf)
         
         bigSummaryDf = rbind(bigSummaryDf, fillerDf)
      }
   }
}

write.csv(bigSummaryDf, file = "bigSummaryDf_BRCA.csv")
save(distributionList, file = "distributionList_BRCA.rda")

library(matrixStats)
library(discover)

cancerList = c("BRCA", "GBM", "LAML", "LGG", "LIHC", "KIRP", "KIRC", "CESC", "STAD", "SKCM", "SARC", "PAAD", "LUSC", "LUAD", "UCEC") 
#cancerList = c("BLCA") 

for(myCancer in cancerList) {
  setwd(paste0("/Users/apple/Desktop/Foundations of Computational Biology and Bioinformatics/FCBB2019Spring/", myCancer, "/"))
  myTab = read.csv(paste0("encoded", myCancer, ".csv")) # change this to appropriate file name  
  #myTab = read.csv(paste0("updated_encoding_", myCancer, ".csv")) # change this to appropriate file name  
    
  rownames(myTab) = myTab$X
  myTab$X = NULL
  
  events <- discover.matrix(myTab)
  
  result.mutex <- pairwise.discover.test(events)
  
  print(result.mutex, fdr.threshold=0.05)
  resultsTable = as.data.frame(result.mutex, q.threshold = 0.05)
  
  file_n = paste0(myCancer, "_mut.csv")
  print(paste0(myCancer, "_mut.csv"))
  write.csv(x = resultsTable, file = file_n)
}

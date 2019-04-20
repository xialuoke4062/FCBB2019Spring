library(discover)
setwd("C:/Users/Dan/Dropbox (CahanLab)/Foundation_FinalProject")

myTab = read.csv("updated_encoding_BLCA.csv")
rownames(myTab) = myTab$X
myTab$X = NULL

events <- discover.matrix(myTab)

result.mutex <- pairwise.discover.test(events)

print(result.mutex, fdr.threshold=0.08)
as.data.frame(result.mutex, q.threshold = 0.08)

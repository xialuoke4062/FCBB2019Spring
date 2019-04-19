library(cancerCellNet)
setwd("/Users/danpeng/Dropbox (CahanLab)/Foundation_FinalProject/brain_tumor/")

oldEncoding = read.csv(file = "LGG_variantClass.csv")

fileName <- "GeneList.txt"
conn <- file(fileName,open="r")
linn <-readLines(conn)
close(conn)

unique_genes = unique(linn)

rownames(oldEncoding) = oldEncoding$X
oldEncoding$X = NULL

oldEncoding = oldEncoding[unique_genes, ]

oldEncoding = na.omit(oldEncoding)


for(colname in colnames(oldEncoding)) {
    tempVector = as.vector(oldEncoding[, colname])
    tempVector[tempVector == "no_mut"] = 0
    tempVector[tempVector == "silent"] = 0
    tempVector[tempVector != 0] = 1
    
    oldEncoding[, colname] = tempVector
}

write.csv(x = oldEncoding, file = "encodedLGG.csv")


## this is to test GBM 
oldEncoding = read.csv(file = "GBM_variantClass.csv")

fileName <- "GeneList.txt"
conn <- file(fileName,open="r")
linn <-readLines(conn)
close(conn)

unique_genes = unique(linn)

rownames(oldEncoding) = oldEncoding$X
oldEncoding$X = NULL

oldEncoding = oldEncoding[unique_genes, ]

oldEncoding = na.omit(oldEncoding)


for(colname in colnames(oldEncoding)) {
    tempVector = as.vector(oldEncoding[, colname])
    tempVector[tempVector == "no_mut"] = 0
    tempVector[tempVector == "silent"] = 0
    tempVector[tempVector != 0] = 1
    
    oldEncoding[, colname] = tempVector
}

write.csv(x = oldEncoding, file = "encodedGBM.csv")
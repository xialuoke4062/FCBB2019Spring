saveFileName = "encodedLUSC.csv"
myMaf = read.csv("/Users/danpeng/Dropbox (CahanLab)/Foundation_FinalProject/GDCdata/GDCdata/TCGA.LUSC.muse.ad676dc8-bc24-440f-bedd-b0885891fb61.DR-10.0.somatic.maf.csv")

genes_row = as.vector(unique(myMaf$Hugo_Symbol))

patient_col = as.vector(unique(myMaf$Tumor_Sample_Barcode))


unique(myMaf$Variant_Classification)

# encode the data using the variation type
tempMatrix = matrix("no_mut", nrow = length(genes_row), ncol = length(patient_col))

rownames(tempMatrix) = genes_row
colnames(tempMatrix) = patient_col

for(patient in patient_col) {
   tempMaf = myMaf[myMaf$Tumor_Sample_Barcode == patient, ]
   tempMatrix[as.vector(tempMaf$Hugo_Symbol), patient] = as.vector(tempMaf$Variant_Classification)
}

oldEncoding = tempMatrix
fileName <- "/Users/danpeng/Desktop/FCBB2019Spring/cleanedData_LGG_GBM/GeneList.txt"
conn <- file(fileName,open="r")
linn <-readLines(conn)
close(conn)

unique_genes = unique(linn)

intersectGenes = intersect(unique_genes, rownames(oldEncoding))
oldEncoding = oldEncoding[intersectGenes, ]

oldEncoding = na.omit(oldEncoding)


for(colname in colnames(oldEncoding)) {
   tempVector = as.vector(oldEncoding[, colname])
   tempVector[tempVector == "no_mut"] = 0
   tempVector[tempVector == "silent"] = 0
   tempVector[tempVector != 0] = 1
   
   oldEncoding[, colname] = tempVector
}

write.csv(x = oldEncoding, file = saveFileName)

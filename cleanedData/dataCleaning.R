myMaf = read.csv("TCGA.GBM.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf.csv")

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

write.csv(tempMatrix, file = "GBM_variantClass.csv")

# to do the LGG group 
myMaf = read.csv("TCGA.LGG.muse.27330980-b14d-4d4d-bd65-5eb19c597c9c.DR-10.0.somatic.maf.csv")

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

write.csv(tempMatrix, file = "LGG_variantClass.csv")
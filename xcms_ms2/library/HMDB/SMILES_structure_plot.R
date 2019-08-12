# BiocManager::install("ChemmineOB")
library(ChemmineOB)
library(ChemmineR)
test_smile = c("CC=O")
test_sdf = smiles2sdf(test_smile)
ChemmineR::plot(test_sdf)

MoNA_MS2_neg = readRDS("../MoNA_MS2_Negative/MoNA_MS2_neg.rds")

spec1 = MoNA_MS2_neg[[1]]
spec2 = MoNA_MS2_neg[[2]]

plot(spec1$spectrum)
spec = spec1
my_SMILES2structure =function(SMILES){
  SDF = smiles2sdf(SMILES)
  ChemmineR::plotStruc(SDF[[1]], regenCoords=T)
}
profvis::profvis({
for(repeating in 1:10){
pdf("testpdf.pdf")
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(spec1$spectrum)
my_SMILES2structure(spec2$SMILES)
plot(spec2$spectrum)
my_SMILES2structure(spec2$SMILES)
dev.off()
}})

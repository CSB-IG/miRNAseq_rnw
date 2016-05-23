library(gProfileR)

enf <- read.delim("matriz_completa_enfermos.txt", header = TRUE, sep = "\t", row.names = 1)

head(enf)
dim(enf)
class(enf)

san <- read.delim("matriz_completa_sanos.txt", header = TRUE, sep = "\t", row.names = 1)

head(san)
dim(san)
class(san)

eset <- cbind(enf, san)
eset <- as.matrix(eset)

head(eset)
dim(eset)
class(eset)

case <- c(1:86)
control <- c(87:172)

length(control)

save(eset, case, control, file = "exp_set.RData")

gconvert <- gconvert(query = rownames(eset), target = "ENSG", mthreshold = 1, filter_na = FALSE)
class(gconvert[,1]) <- "integer"
gconvert <- gconvert[sort.list(gconvert[,1]),]
genesymbol <- gconvert[,4]
genesymbol <- ifelse(genesymbol == "N/A", as.vector(rownames(eset)), as.vector(genesymbol))

TF <- read.table(file= "../annotation/TF_vaquerizas.txt", )
TF <- as.character(TF[[1]])

length(TF)

TFinGS <- which(genesymbol %in% TF)
mirnas <- grep("hsa-", genesymbol)

regulators <- rownames(eset)[c(mirnas, TFinGS)]

head(regulators)

write.table(regulators, file = "regulators.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
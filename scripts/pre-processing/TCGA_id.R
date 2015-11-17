RNAseqV2 <- read.csv("rnaseq.csv")
id <- RNAseqV2[,1]
anyDuplicated(id)
id4 <- gsub('.{4}$', '', id)
duplicated(id4)
sum(duplicated(id4))

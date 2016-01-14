

library(AIMS)
# Arguments
# eset = matrix with genes as rows and samples columns
# EntrezID = list of length same as rownames eset with gene Entrez ID in order

dataset <- read.table("RNAseq_752_norm.rsem.aracne.txt", header=T, sep="\t", row.names=1)
dataset <- data.matrix(dataset)
dannot <- read.table("annot.txt", header=F, sep="\t", row.names=1)
colnames(dannot)<-c("probe", "EntrezGene.ID")

genes <- rownames(dataset)
IDs <- dannot[genes, ]
genes <- IDs[, 2]

Preds <- applyAIMS(eset=dataset, EntrezID=genes)

names(Preds)

table(PAM50Preds$cl)


write.table(PAM50Preds$cl, quote = F, sep = "\t", row.names = T, file = "AIMS_RNAseq.txt")

# table(PAM50Preds$cl)

# Basal   Her2   LumA   LumB Normal 
#   131     75    288    198     60 
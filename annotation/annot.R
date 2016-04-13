library(gProfileR)

dannot <- read.table("annot.txt", header=F, sep="\t", row.names=1)
RNA86 <- scan(file="probes_RNAseq_86", what="1", sep="\t")
RNA752 <- scan(file="probes_RNAseq_752", what="1", sep="\t")
annot86 <- dannot[RNA86,]
annot752 <- dannot[RNA752,]
a86 <- annot86[,2]
a752 <- annot752[,2]
ensg86 <- gconvert(query=a86, numeric_ns="ENTREZGENE_ACC", target="ENSG", mthreshold=1, filter_na=FALSE)
ensg752 <- gconvert(query=a752, numeric_ns="ENTREZGENE_ACC", target="ENSG", mthreshold=1, filter_na=FALSE)
rownames(ensg86)<- ensg86[,1]
rownames(ensg752)<- ensg752[,1]
ensg86 <- ensg86[order(as.numeric(as.character(rownames(ensg86)))), ]
ensg752 <- ensg752[order(as.numeric(as.character(rownames(ensg752)))), ]
annot86 <- cbind(annot86, ensg86[, "target"])
annot752 <- cbind(annot752, ensg752[, "target"])
vaq <- scan(file="TF.txt", what="1", sep="\t")
annot86a <- subset(annot86, annot86[,3] %in% vaq)
annot752a <- subset(annot752, annot752[,3] %in% vaq)
annot86a <- annot86a[sort(rownames(annot86a)), ]
annot752a <- annot752a[sort(rownames(annot752a)), ]
write.table(annot86a, quote = F, sep = "\t", row.names = T, col.names = NA, file = "annot_86_RNAseq_TF.txt")
write.table(annot752a, quote = F, sep = "\t", row.names = T, col.names = NA, file = "annot_752_RNAseq_TF.txt")
TF86 <- cbind(intersect(RNA86, rownames(annot86a)), 2)
TF572 <- cbind(intersect(RNA752, rownames(annot752a)), 2)
genes86 <- cbind(setdiff(RNA86, rownames(annot86a)), 1)
genes752 <- cbind(setdiff(RNA752, rownames(annot752a)), 1)
write.table(TF86, quote = F, sep = "\t", row.names = F, col.names = F, file = "RNAseq_86_TF.txt")
write.table(TF752, quote = F, sep = "\t", row.names = F, col.names = F, file = "RNAseq_752_TF.txt")
write.table(genes86, quote = F, sep = "\t", row.names = F, col.names = F, file = "RNAseq_86_genes.txt")
write.table(genes752, quote = F, sep = "\t", row.names = F, col.names = F, file = "RNAseq_752_genes.txt")

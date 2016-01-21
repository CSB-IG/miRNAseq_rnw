

genes <- read.table("RNAseq_752_norm.rsem.aracne.txt", sep="\t", header=T, row.names=1)

miRNAs <- read.table("miRNA_752_normTMM.txt", sep="\t", header=T, row.names=1)

> a <- colnames(genes)
> b <- colnames(miRNAs)
> identical(a, b)
[1] TRUE

> red <- rbind(miRNAs, genes)
dim(red)

write.table(red, quote = F, sep = "\t", row.names = T, col.names = NA, file = "matriz_completa.txt")
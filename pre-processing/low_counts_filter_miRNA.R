# x <- expression matrix mature miRNAseq counts

# sum values per row before filtering
# sums <- rowSums(x)

# filter targets with a minimum of  5  counts  in  at  least  25%  of  the  samples(752) = 188
miRNA <- x[rowSums(x) >= 188, , drop=FALSE]

# > dim(miRNA)
# [1] 715 752

# sums for the remaining genes
# f <- rowSums(miRNA)

# save final expression matrix
write.table(miRNA, quote = F, sep = "\t", row.names = T, col.names = NA, file = "miRNAs_752_maduros_filtro.txt")

# rowSums simple example:

#mymat <- matrix(sample(c(1,0), 100, r=TRUE), ncol=10)
#rowSums(mymat)
#mymat[rowSums(mymat) == 5, , drop=FALSE]
#s <- mymat[rowSums(mymat) >= 5, , drop=FALSE]
#rowSums(s)

# matrix format: rownames: GENE/miR ID, colnames: patient code: XXXX.01(case)/11(control)

# raw counts 
x <- read.table(file="expn_matrix_mimat.txt", sep="\t", header=T, row.names=1)
> dim(x)
[1] 2588  753
names <- colnames(x)
names1 <- gsub('.{13}$', '', names)
names2 <- substring(names1, 9)
colnames(x) <- names2
x$A245.01 <- NULL # only for 753 samples

# save final expression matrix
write.table(x, quote = F, sep = "\t", row.names = T, col.names = NA, file = "miRNAs_752_maduros.txt")

# normalized counts
xn <- read.table(file="expn_matrix_mimat_norm.txt", sep="\t", header=T, row.names=1)
> dim(xn)
[1] 2588  753
nnames <- colnames(xn)
nnames1 <- gsub('.{13}$', '', nnames)
nnames2 <- substring(nnames1, 9)
colnames(xn) <- nnames2


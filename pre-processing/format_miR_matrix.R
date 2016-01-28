# matrix format: rownames: GENE/miR ID, colnames: patient code: XXXX.01(case)/11(control)

# raw counts 
x <- read.table(file="expn_matrix_mimat.txt", sep="\t", header=T, row.names=1)
> dim(x)
[1] 2588  753
names <- colnames(x)
names1 <- gsub('.{13}$', '', names)
names2 <- substring(names1, 9)
colnames(x) <- names2
x$A245.01 <- NULL

# save final expression matrix
write.table(x, quote = F, sep = "\t", row.names = T, col.names = NA, file = "miRNAs_752_maduros.txt")



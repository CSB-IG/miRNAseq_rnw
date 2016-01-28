# x <- expression matrix mature miRNAseq counts
# enfermos <- expression matrix mature miRNAseq counts for a group (cases)
# sanos <- expression matrix mature miRNAseq counts for a group (controls)

# sum values per row before filtering
# sums <- rowSums(x)

# filter targets with a minimum of  5  counts  in  at  least  25%  of  the  samples(752) = 188
miRNA <- x[rowSums(x) >= 188, , drop=FALSE]

# filter targets with a minimum of  5  counts  in  at  least  25%  of  the  samples(86) = 22
enfermos <- enfermos[rowSums(enfermos) >= 22, , drop=FALSE]

# filter targets with a minimum of  5  counts  in  at  least  25%  of  the  samples(86) = 22
sanos <- sanos[rowSums(sanos) >= 22, , drop=FALSE]

# > dim(enfermos)                                                                                                                       
# [1] 690  86
# > dim(sanos)
# [1] 658  86
# > length(intersect(rownames(enfermos), rownames(sanos)))
# [1] 633

# set two groups to smaller gene number(658):
enfermos <- enfermos[intersect(rownames(enfermos), rownames(sanos)), ]
sanos <- sanos[intersect(rownames(enfermos), rownames(sanos)), ]
identical(rownames(sanos), rownames(enfermos))
identical(gsub('.{3}$', '',colnames(enfermos)), gsub('.{3}$', '',colnames(sanos)))

# save final expression matrix
write.table(enfermos, quote = F, sep = "\t", row.names = T, col.names = NA, file = "miRNAs_enfermos_maduros_filtro.txt")
write.table(sanos, quote = F, sep = "\t", row.names = T, col.names = NA, file = "miRNAs_sanos_maduros_filtro.txt")

# rowSums simple example:

#mymat <- matrix(sample(c(1,0), 100, r=TRUE), ncol=10)
#rowSums(mymat)
#mymat[rowSums(mymat) == 5, , drop=FALSE]
#s <- mymat[rowSums(mymat) >= 5, , drop=FALSE]
#rowSums(s)

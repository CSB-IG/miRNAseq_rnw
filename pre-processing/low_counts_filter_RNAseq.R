# read expression matrix
RNAseq <- read.csv("RNAseq_norm_casos_752.csv")

# sum values per row before filtering
# sums <- rowSums(x = RNAseq)

# filtering sum >= 7520 / 10 TPM
RNAseq <- RNAseq[rowSums(x = RNAseq) >= 7520, , drop=FALSE]

# sums for the remaining genes
# f <- rowSums(x = filtro)

# GENEID index corresponding to remaining row.names
indexID <- as.numeric(row.names(RNAseq))

# colnames(indexID)[1] = "GENEindex"

# complete GENEID gene names list
ID <- read.table("GeneID", header = TRUE, stringsAsFactors = FALSE)

# filter GENEID list by GENEindex
GENEID_filtro <- data.frame(ID[indexID,])

# save and ignore first 10 genes from GENEID
write.table(GENEID_filtro[11:15982,], row.names = F, col.names = F, quote = F, sep = "\t", file = "GENEID_filtro.txt")
# edit in emacs
# non-unique value when setting 'row.names': ‘SLC35E2'
# SLC35E2|728661 -> SLC35E2B
# SLC35E2|9906   -> SLC35E2


# filtered expression matrix with GENEIDs, ignoring 10 genes without GENEID
m <- RNAseq[11:15752,]
ID <- scan("GENEID_filtro.txt", what = "", sep="\t")
row.names(m) <- ID

# save final expression matrix
write.table(filtro[11:15982,], quote = F, sep = "\t", row.names = F, col.names = T, file = "RNAseq_752_filtered.txt")

# other options 
# to create a dataframe from indexes
# indexID <- data.frame(as.numeric(row.names(RNAseq)))

# make filter's first column the ID in data.frame
# indexID <- data.frame(as.numeric(row.names(filtro)), filtro)

# rowSums simple example:

#mymat <- matrix(sample(c(1,0), 100, r=TRUE), ncol=10)
#rowSums(mymat)
#mymat[rowSums(mymat) == 5, , drop=FALSE]
#s <- mymat[rowSums(mymat) >= 5, , drop=FALSE]
#rowSums(s)

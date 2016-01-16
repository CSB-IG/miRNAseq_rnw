# read expresssion matrix
RNAseq <- read.csv("RNAseq_norm_casos_86.csv")

# sum values per row before filtering
# sums <- rowSums(x = RNAseq)

# filtering suma >= 7520/ 10 TPM
RNAseq <- RNAseq[rowSums(x = RNAseq) >= 860, , drop=FALSE]

# sums for the remaining genes 
# f <- rowSums(x = RNAseq)

# GENEID index corresponding to remaining row.names
indexID <- as.numeric(row.names(RNAseq))

# colnames(indexID)[1] = "GENEindex"

# complete GENEID gene names list 
ID <- read.table("GeneID", header = TRUE, stringsAsFactors = FALSE)

# filter GENEID list by GENEindex 
GENEID_filtro <- data.frame(ID[indexID,])

# save and ignore first 10 genes from GENEID
write.table(GENEID_filtro[11:15752,], row.names = F, col.names = F, quote = F, sep = "\t", file = "GENEID_enfermos.txt")
# edit in emacs:
# non-unique value when setting 'row.names': ‘SLC35E2'
# SLC35E2|728661 -> SLC35E2B
# SLC35E2|9906   -> SLC35E2

# filtered expression matrix with GENEIDs, ignoring first 10 genes without GENEID
m <- RNAseq[11:15752,] 
ID <- scan("GENEID_filtro.txt", what = "", sep="\t")
row.names(m) <- ID

# save final expression matrix
write.table(m, quote = F, sep = "\t", row.names = T, col.names = NA, file = "enfermos.86.txt")

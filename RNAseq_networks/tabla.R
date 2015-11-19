
# read expression matrix
x <- read.table("RNAseq_752_filtered.txt", header = T)

# gene ID list
# the gene ID: SLC35E2 entry is repeated, modify first SLC35E2|728661 entry -> to SLC35E2B
# for GENEID names (if header): ID <- scan("GENEID_filtro.txt", what = "", sep="\t", skip = 1)
ID <- scan("GENEID_filtro.txt", what = "", sep="\t")

# add ID as rownames
row.names(x) <- (ID)

# save expression matrix with corresponding GENEIDs
write.table(x, file = "RNAseq_752.15972.rsem.me.txt", sep = "\t", row.names = T, quote = F, col.names = NA)

# compare gene ID RNAseq_752.15972.rsem.me.txt vs RNAseq_752_norm.rsem.aracne.txt
# library(gplots)
# y <- read.table("RNAseq_752_norm.rsem.aracne.txt", header = T, row.names = 1)
# a <- row.names(y)
# b <- ID
# c <- venn(list(hugo = a,  diana = b), show.plot = T)


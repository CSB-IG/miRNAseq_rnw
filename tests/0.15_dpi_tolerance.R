
# load p-value prunned network matrices without any dpi prunning

sanos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/sanos_pvalue_full_adjmtx.txt"))
enfermos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/enfermos_pvalue_full_adjmtx.txt"))

colnames(sanos) <- rownames(sanos)
colnames(enfermos) <- rownames(enfermos)

dim(sanos)
dim(enfermos)

isSymmetric(sanos)
isSymmetric(enfermos)

max(sanos)
max(enfermos)

library(minet)

sanos_dpi <- aracne(mim = sanos, eps = 0.15)
enfermos_dpi <- aracne(mim = enfermos, eps = 0.15)

max(sanos_dpi)
max(enfermos_dpi)

write.table(sanos_dpi, file = "sanos_pvalue_full_adjmtx_dpi_015.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
write.table(enfermos_dpi, file = "enfermos_pvalue_full_adjmtx_dpi_015.txt", sep = "\t", col.names= NA, row.names= T, quote = F)


redes <- function(g){
	library(igraph)
	g1 <- graph.adjacency(adjmatrix= g, mode='undirected', diag=F, weighted=T)
	g2 <- get.data.frame(g1)
	return(g2)
}

red_sanos <- redes(sanos_dpi)
red_enfermos <- redes(enfermos_dpi)


write.table(red_sanos, file = "sanos_pvalue_full_sif_dpi_015.txt", sep = "\t", col.names= T, row.names= F, quote = F )
write.table(red_enfermos, file = "enfermos_pvalueq_sif_full_dpi_015.txt", sep = "\t", col.names= T, row.names= F, quote = F )



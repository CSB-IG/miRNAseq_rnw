
library(igraph)

# sacar dpi a p-value 1
 
sanos <- read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/sanos_p1_full_norm_adjmtx.txt")
enfermos <- read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/enfermos_p1_full_norm_adjmtx.txt")

sanos <- data.matrix(sanos)
enfermos <- data.matrix(enfermos)

colnames(sanos) <- rownames(sanos)
colnames(enfermos) <- rownames(enfermos)

dim(sanos)
dim(enfermos)

isSymmetric(sanos)
isSymmetric(enfermos)

max(sanos)
max(enfermos)

library(minet)
sanos_dpi <- aracne(mim = sanos, eps = 0.1)
enfermos_dpi <- aracne(mim = enfermos, eps = 0.1)

write.table(sanos_dpi, file = "sanos_p1_full_adjmtx_dpi.txt", sep = "\t", col.names= NA, row.names= T, quote = F)

sg <- graph.adjacency(adjmatrix= sanos_dpi, mode='undirected', diag=F, weighted=T)
sifgs <- get.data.frame(sg)
write.table(sifgs, file = "sanos_p1_full_sif_dpi.txt", sep = "\t", col.names= T, row.names= F, quote = F )

write.table(enfermos_dpi, file = "enfermos_p1_full_adjmtx_dpi.txt", sep = "\t", col.names= NA, row.names= T, quote = F)

eg <- graph.adjacency(adjmatrix= enfermos_dpi, mode='undirected', diag=F, weighted=T)
sifge <- get.data.frame(eg)
write.table(sifge, file = "enfermos_p1_sif_full_dpi.txt", sep = "\t", col.names= T, row.names= F, quote = F )

max(sanos_dpi)
max(enfermos_dpi)

isSymmetric(sanos_dpi)
isSymmetric(enfermos_dpi)

# calcular p value y corte

mir <- colnames(enfermos)[grepl("hsa*", colnames(enfermos))]
head(mir)
length(mir)

genes <- setdiff(colnames(enfermos), mir)
head(genes)
length(genes)

sanos_q <- sanos[upper.tri(sanos, diag=FALSE)]
quantile(sanos_q)

enfermos_q <- enfermos[upper.tri(enfermos, diag=FALSE)]
quantile(enfermos_q)


# filtrar por p-value




# guardar

write.table(mir_gen, file = "sanos_mir_gen_p1_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xmg <- graph.adjacency(adjmatrix= mir_gen, mode='undirected', diag=F, weighted=T)
xsifmg <- get.data.frame(xmg)
write.table(xsifmg, file = "sanos_mir_gen_p1_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )
 

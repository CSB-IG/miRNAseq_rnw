# sin ordenar x
# solo agregar nombres de columnas

x[is.na(x)] <- 0
class(x) <- "numeric"
names <- gsub('.{6}$','', adjs)
colnames(x) <- names

a <- x[15593:16225, 15593:16225]
> dim(a)
[1] 633 633
> sum(is.na(a))
[1] 400689
> 633**2
[1] 400689

# b no tiene miRNAS en rownames

b <- x[1:15592, 1:15592]

adjs_e <- list.files(path="~/TCGA_miRNA_BC/miRNAs_resultados/miRNAS_86_resultados/RNAseq_casos_resultados_ARACNE", pattern="*_1.adj")
length(adjs_e)
# [1] 15742
length(intersect(adjs, adjs_e))
# [1] 15136
genes <- (intersect(adjs, adjs_e))
genes <- gsub('.{6}$','', genes)
b <- b[genes, ]

b <- b[sort(rownames(b)), sort(colnames(b)) ]
length(setdiff(rownames(b), colnames(b)))
length(setdiff(colnames(b), rownames(b)))
# ahora genes tiene los miRNAS y los genes que hay que quitar
dif <- setdiff(colnames(b), rownames(b))
mir <- dif[grepl("hsa*", dif)]
y <- t(b[,mir])
y2 <- x[mir,mir]

identical(colnames(y2), rownames(y))
identical(rownames(y2), rownames(y))
y <- cbind(y, y2)
faltan <- setdiff(colnames(y), colnames(b))

m <- x[genes,faltan]
genes_faltan <- (setdiff(genes, faltan))
final <- c(genes_faltan, mir)
m <- m[sort(rownames(m)),]
b <- b[sort(rownames(b)), final]
identical(rownames(b), rownames(m))
b <- cbind(b,m)
b <- b[,sort(colnames(b))]
y <- y[,sort(colnames(y))]
identical(colnames(b), colnames(y))
enfermos <-rbind(b,y)
enfermos <- enfermos[sort(rownames(enfermos)), sort(colnames(enfermos)) ]

normalit<-function(m){
     (m - min(m))/(max(m)-min(m))
}

normalit(m)

write.table(sanos_norm, file = "sanos_p1_full_norm_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xg_norm <- graph.adjacency(adjmatrix= sanos_norm, mode='undirected', diag=F, weighted=T)
xsif_norm <- get.data.frame(xg_norm)
write.table(xsif_norm, file = "sanos_p1_full_norm_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

mirnas <- sanos[mir, mir]
write.table(mirnas, file = "sanos_mir_mir_p1_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xgmm <- graph.adjacency(adjmatrix= mirnas, mode='undirected', diag=F, weighted=T)
xsifmm <- get.data.frame(xgmm)
write.table(xsifmm, file = "sanos_mir_mir_p1_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

mirnas_norm <- sanos_norm[mir, mir]
write.table(mirnas_norm, file = "sanos_mir_mir_p1_norm_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xgmm_norm <- graph.adjacency(adjmatrix= mirnas_norm, mode='undirected', diag=F, weighted=T)
xsifmm_norm <- get.data.frame(xgmm_norm)
write.table(xsifmm_norm, file = "sanos_mir_mir_p1_norm_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

gengen <- sanos[genes, genes]
write.table(gengen, file = "sanos_gen_gen_p1_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xgg <- graph.adjacency(adjmatrix= gengen, mode='undirected', diag=F, weighted=T)
xsifgg <- get.data.frame(xgg)
write.table(xsifgg, file = "sanos_gen_gen_p1_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

gengen_norm <- sanos_norm[genes, genes]
write.table(gengen_norm, file = "sanos_gen_gen_p1_norm_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xgg_norm <- graph.adjacency(adjmatrix= gengen_norm, mode='undirected', diag=F, weighted=T)
xsifgg_norm <- get.data.frame(xgg_norm)
write.table(xsifgg_norm, file = "sanos_gen_gen_p1_norm_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

mir_gen <- sanos
mir_gen[genes,genes] <- NA
mir_gen[mir,mir] <- NA
mir_gen[is.na(mir_gen)] <- 0
isSymmetric(mir_gen)

write.table(mir_gen, file = "sanos_mir_gen_p1_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xmg <- graph.adjacency(adjmatrix= mir_gen, mode='undirected', diag=F, weighted=T)
xsifmg <- get.data.frame(xmg)
write.table(xsifmg, file = "sanos_mir_gen_p1_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

mir_gen_norm <- sanos_norm
mir_gen_norm[genes,genes] <- NA
mir_gen_norm[mir,mir] <- NA
mir_gen_norm[is.na(mir_gen)] <- 0
isSymmetric(mir_gen_norm)

write.table(mir_gen_norm, file = "sanos_mir_gen_p1_norm_adjmtx.txt", sep = "\t", col.names= NA, row.names= T, quote = F)
xmg_norm <- graph.adjacency(adjmatrix= mir_gen_norm, mode='undirected', diag=F, weighted=T)
xsifmg_norm <- get.data.frame(xmg_norm)
write.table(xsifmg_norm, file = "sanos_mir_gen_p1_sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

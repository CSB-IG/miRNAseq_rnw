
enfermos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/enfermos_p1_full_adjmtx.txt"))
colnames(enfermos) <- rownames(enfermos)
dim(enfermos)
isSymmetric(enfermos)

sanos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/sanos_p1_full_adjmtx.txt"))
colnames(sanos) <- rownames(sanos)
dim(sanos)
isSymmetric(sanos) 

mir <- rownames(enfermos)[grepl("hsa*", rownames(enfermos))]
genes <- rownames(enfermos)[grepl("^(?!hsa*)", rownames(enfermos), perl=T)]

num_int_genes <- ((length(genes)**2)-length(genes))/2
num_int_mir <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 

(100/num_int_genes)*25334
100-((100/num_int_genes)*25334)

(100/num_int_mir)*14892
100-((100/num_int_mir)*14892)

# threshold for genes from miRNA int number

# > (100/num_int_genes)*25334
# [1] 0.02211771
# > 100-((100/num_int_genes)*25334)                                                                                                                                                                                                                                                                                              
# [1] 99.97788


# threshold for miRNAs from gene int number

# > num_int_mir <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 
# > (100/num_int_mir)*14892
# [1] 0.1522526
# > 100-((100/num_int_mir)*14892)
# [1] 99.84775

# quantiles for 2533 miR and 1489 gene int

(100/num_int_mir)*2533
100-((100/num_int_mir)*2533)

(100/num_int_genes)*1489
100-((100/num_int_genes)*1489)

# > (100/num_int_mir)*2533
# [1] 0.02589684
# > 100-((100/num_int_mir)*2533)
# [1] 99.9741
# > 
# > (100/num_int_genes)*1489
# [1] 0.001299963
# > 100-((100/num_int_genes)*1489)
# [1] 99.9987

# quantiles for 253340 miR and 148920 gene int

(100/num_int_mir)*253340
100-((100/num_int_mir)*253340)

(100/num_int_genes)*148920
100-((100/num_int_genes)*148920)

# > (100/num_int_mir)*253340
# [1] 2.590093
# > 100-((100/num_int_mir)*253340)
# [1] 97.40991
# > 
# > (100/num_int_genes)*148920
# [1] 0.1300138
# > 100-((100/num_int_genes)*148920)
# [1] 99.86999

# thresholding

mi_mirnas <- function(x){
	return(c((x[mir,mir][upper.tri(x[mir, mir])]), as.vector(x[mir, genes])))
}

mi_genes <- function(x){
	return(x[genes,genes][upper.tri(x[genes, genes])])
}

quantile2mi <- function(q,v){
	return(quantile(v, probs=q, na.rm=T))
}

matriz_genes <- function(x,n){
	x = x[genes,genes]
	x[x < n] <- 0
	return(x)
}


matriz_mir <- function(x,n){
	x[genes,genes] <- 0
	x[x < n] <- 0
	return(x)
}

enfermos_mir_q <- mi_mirnas(enfermos)
enfermos_genes_q <- mi_genes(enfermos)

sanos_mir_q <- mi_mirnas(sanos)
sanos_genes_q <- mi_genes(sanos)

# mismo numero int en genes que int miRNAs
enfermos_thresholdmir_mir <- matriz_mir(enfermos, (as.vector(quantile2mi(0.99741, enfermos_mir_q))))
enfermos_thresholdmir_gen <- matriz_genes(enfermos, (as.vector(quantile2mi(0.99978, enfermos_genes_q))))

enfermos_thresholdmir_mir[genes,genes] <- enfermos_thresholdmir_gen

sanos_thresholdmir_mir <- matriz_mir(sanos, (as.vector(quantile2mi(0.99741, sanos_mir_q))))
sanos_thresholdmir_gen <- matriz_genes(sanos, (as.vector(quantile2mi(0.99978, sanos_genes_q))))

sanos_thresholdmir_mir[genes,genes] <- sanos_thresholdmir_gen

# mismo numero int en miRNAS que int genes
enfermos_thresholdgen_mir<- matriz_mir(enfermos, (as.vector(quantile2mi(0.99848, enfermos_mir_q))))
enfermos_thresholdgen_gen <- matriz_genes(enfermos, (as.vector(quantile2mi(0.99987, enfermos_genes_q))))

enfermos_thresholdgen_mir[genes,genes] <- enfermos_thresholdgen_gen

sanos_thresholdgen_mir <- matriz_mir(sanos, (as.vector(quantile2mi(0.99848, sanos_mir_q))))
sanos_thresholdgen_gen <- matriz_genes(sanos, (as.vector(quantile2mi(0.99987, sanos_genes_q))))

sanos_thresholdgen_mir[genes,genes] <- sanos_thresholdgen_gen

# orden de mag menor
enfermos_thresholddown_mir <- matriz_mir(enfermos, (as.vector(quantile2mi(0.999741, enfermos_mir_q))))
enfermos_thresholddown_gen <- matriz_genes(enfermos, (as.vector(quantile2mi(0.999987, enfermos_genes_q))))

enfermos_thresholddown_mir[genes,genes] <- enfermos_thresholddown_gen

sanos_thresholddown_mir <- matriz_mir(sanos, (as.vector(quantile2mi(0.999741, sanos_mir_q))))
sanos_thresholddown_gen <- matriz_genes(sanos, (as.vector(quantile2mi(0.999987, sanos_genes_q))))

sanos_thresholddown_mir[genes,genes] <- sanos_thresholddown_gen

# orden de mag mayor
enfermos_thresholdup_mir <- matriz_mir(enfermos, (as.vector(quantile2mi(0.9741, enfermos_mir_q))))
enfermos_thresholdup_gen <- matriz_genes(enfermos, (as.vector(quantile2mi(0.99870, enfermos_genes_q))))

enfermos_thresholdup_mir[genes,genes] <- enfermos_thresholdup_gen

sanos_thresholdup_mir <- matriz_mir(sanos, (as.vector(quantile2mi(0.9741, sanos_mir_q))))
sanos_thresholdup_gen <- matriz_genes(sanos, (as.vector(quantile2mi(0.99870, sanos_genes_q))))

sanos_thresholdup_mir[genes,genes] <- sanos_thresholdup_gen

# guardar

write.table(enfermos_thresholdmir_mir, file="enfermos_thresholdmir.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholdmir_mir, file="sanos_thresholdmir.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(enfermos_thresholdgen_mir, file="enfermos_thresholdgen.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholdgen_mir, file="sanos_thresholdgen.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(enfermos_thresholddown_mir, file="enfermos_thresholddown.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholddown_mir, file="sanos_thresholddown.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(enfermos_thresholdup_mir, file="enfermos_thresholdup.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholdup_mir, file="sanos_thresholdup.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)


# dpi
library(minet)

enfermos_thresholdmir_minet <- aracne(mim = enfermos_thresholdmir_mir, eps = 0.1)
sanos_thresholdmir_minet <- aracne(mim = sanos_thresholdmir_mir, eps = 0.1)

enfermos_thresholdgen_minet <- aracne(mim = enfermos_thresholdgen_mir, eps = 0.1)
sanos_thresholdgen_minet <- aracne(mim = sanos_thresholdgen_mir, eps = 0.1)

enfermos_thresholddown_minet <- aracne(mim = enfermos_thresholddown_mir, eps = 0.1)
sanos_thresholddown_minet <- aracne(mim = sanos_thresholddown_mir, eps = 0.1)

enfermos_thresholdup_minet <- aracne(mim = enfermos_thresholdup_mir, eps = 0.1)
sanos_thresholdup_minet <- aracne(mim = sanos_thresholdup_mir, eps = 0.1)

# save

write.table(enfermos_thresholdmir_minet, file="enfermos_thresholdmir_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholdmir_minet, file="sanos_thresholdmir_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(enfermos_thresholdgen_minet, file="enfermos_thresholdgen_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholdgen_minet, file="sanos_thresholdgen_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(enfermos_thresholddown_minet, file="enfermos_thresholddown_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholddown_minet, file="sanos_thresholddown_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(enfermos_thresholdup_minet, file="enfermos_thresholdup_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_thresholdup_minet, file="sanos_thresholdup_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

isSymmetric(enfermos_thresholdmir_minet)
isSymmetric(sanos_thresholdmir_minet)
isSymmetric(enfermos_thresholdgen_minet)
isSymmetric(sanos_thresholdgen_minet)
isSymmetric(enfermos_thresholddown_minet)
isSymmetric(sanos_thresholddown_minet)
isSymmetric(enfermos_thresholdup_minet)
isSymmetric(sanos_thresholdup_minet)

redes <- function(g){
	library(igraph)
	g1 <- graph.adjacency(adjmatrix= g, mode='undirected', diag=F, weighted=T)
	g2 <- get.data.frame(g1)
	return(g2)
}

red_enfermos_thresholdmir <- redes(enfermos_thresholdmir_minet)
red_sanos_thresholdmir <- redes(sanos_thresholdmir_minet)

red_enfermos_thresholdgen <- redes(enfermos_thresholdgen_minet)
red_sanos_thresholdgen <- redes(sanos_thresholdgen_minet)

red_enfermos_thresholddown <- redes(enfermos_thresholddown_minet)
red_sanos_thresholddown <- redes(sanos_thresholddown_minet)

red_enfermos_thresholdup <- redes(enfermos_thresholdup_minet)
red_sanos_thresholdup <- redes(sanos_thresholdup_minet)

# save

write.table(red_enfermos_thresholdmir, file="red_enfermos_thresholdmir_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(red_sanos_thresholdmir, file="red_sanos_thresholdmir_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(red_enfermos_thresholdgen, file="red_enfermos_thresholdgen_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(red_sanos_thresholdgen, file="red_sanos_thresholdgen_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(red_enfermos_thresholddown, file="red_enfermos_thresholddown_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(red_sanos_thresholddown, file="red_sanos_thresholddown_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

write.table(red_enfermos_thresholdup, file="red_enfermos_thresholdup_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(red_sanos_thresholdup, file="red_sanos_thresholdup_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)


# > dim(red_enfermos_thresholdmir)
# [1] 43888     3
# > dim(red_sanos_thresholdmir)
# [1] 48060     3
# > dim(red_enfermos_thresholdgen)
# [1] 27031     3
# > dim(red_sanos_thresholdgen)
# [1] 28907     3
# > dim(red_enfermos_thresholddown)
# [1] 3877    3
# > dim(red_sanos_thresholddown)
# [1] 4018    3
# > dim(red_enfermos_thresholdup)
# [1] 321099      3
# > dim(red_sanos_thresholdup)
# [1] 328648      3

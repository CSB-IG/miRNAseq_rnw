enfermos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/enfermos_p1_full_norm_adjmtx.txt"))
colnames(enfermos) <- rownames(enfermos)
dim(enfermos)
isSymmetric(enfermos)

sanos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/sanos_p1_full_norm_adjmtx.txt"))
colnames(sanos) <- rownames(sanos)
dim(sanos)
isSymmetric(sanos) 

enfermos_sn <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/enfermos_p1_full_adjmtx.txt"))
colnames(enfermos_sn) <- rownames(enfermos_sn)
dim(enfermos_sn)
isSymmetric(enfermos_sn)

sanos_sn <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/sanos_p1_full_adjmtx.txt"))
colnames(sanos_sn) <- rownames(sanos_sn)
dim(sanos_sn)
isSymmetric(sanos_sn)

> identical(sanos, sanos_sn)
[1] FALSE
> identical(colnames(sanos), colnames(sanos_sn))                                                                                                                          
[1] TRUE
> identical(rownames(sanos), rownames(sanos_sn))
[1] TRUE

scale_mtx <-function(m){
	if(min(m[upper.tri(m)]) == 0){
		a <- (m-min(m))/(max(m)-min(m))
		return(a)
	}else{
		a <- (m - min(m[upper.tri(m)]))/(max(m)-min(m[upper.tri(m)]))
		diag(a) <- 0
		return(a)
		}
}

nuevo_sanos <- scale_matrix(sanos_sn)
nuevo_enfermos <- scale_matrix(enfermos_sn)


nuevo_sanos_round <- round(nuevo_sanos, digits = 6)
nuevo_enfermos_round <- round(nuevo_enfermos, digits = 6)

sanos_round <- round(sanos, digits = 6)
enfermos_round <- round(enfermos, digits = 6)

identical(nuevo_sanos_round,sanos_round)
identical(nuevo_enfermos_round,enfermos_round)
[1] TRUE


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

mir <- rownames(enfermos)[grepl("hsa*", rownames(enfermos))]
genes <- rownames(enfermos)[grepl("^(?!hsa*)", rownames(enfermos), perl=T)]



enfermos_mir_q <- mi_mirnas(enfermos)
enfermos_genes_q <- mi_genes(enfermos)

sanos_mir_q <- mi_mirnas(sanos)
sanos_genes_q <- mi_genes(sanos)




enfermos_mir_q_nuevo <- mi_mirnas(nuevo_enfermos)
enfermos_genes_q_nuevo <- mi_genes(nuevo_enfermos)

sanos_mir_q_nuevo <- mi_mirnas(nuevo_sanos)
sanos_genes_q_nuevo <- mi_genes(nuevo_sanos)




enfermos_mir_q_round <- mi_mirnas(enfermos_round)
enfermos_genes_q_round <- mi_genes(enfermos_round)

sanos_mir_q_round <- mi_mirnas(sanos_round)
sanos_genes_q_round <- mi_genes(sanos_round)



enfermos_mir_q_nuevo_round <- mi_mirnas(nuevo_enfermos_round)
enfermos_genes_q_nuevo_round <- mi_genes(nuevo_enfermos_round)

sanos_mir_q_nuevo_round <- mi_mirnas(nuevo_sanos_round)
sanos_genes_q_nuevo_round <- mi_genes(nuevo_sanos_round)


# matriz anterior
enfermos_mir_viejo <- matriz_mir(enfermos, (as.vector(quantile2mi(0.99741, enfermos_mir_q))))
enfermos_gen_viejo <- matriz_genes(enfermos, (as.vector(quantile2mi(0.99987, enfermos_genes_q))))

enfermos_mir_viejo[genes,genes] <- enfermos_gen_viejo

sanos_mir_viejo <- matriz_mir(sanos, (as.vector(quantile2mi(0.99741, sanos_mir_q))))
sanos_gen_viejo <- matriz_genes(sanos, (as.vector(quantile2mi(0.99987, sanos_genes_q))))

sanos_mir_viejo[genes,genes] <- sanos_gen_viejo

# matriz nueva
enfermos_mir_nuevo <- matriz_mir(nuevo_enfermos, (as.vector(quantile2mi(0.99741, enfermos_mir_q_nuevo))))
enfermos_gen_nuevo <- matriz_genes(nuevo_enfermos, (as.vector(quantile2mi(0.99987, enfermos_genes_q_nuevo))))

enfermos_mir_nuevo[genes,genes] <- enfermos_gen_nuevo

sanos_mir_nuevo <- matriz_mir(nuevo_sanos, (as.vector(quantile2mi(0.99741, sanos_mir_q_nuevo))))
sanos_gen_nuevo <- matriz_genes(nuevo_sanos, (as.vector(quantile2mi(0.99987, sanos_genes_q_nuevo))))

sanos_mir_nuevo[genes,genes] <- sanos_gen_nuevo



# viejo round 
enfermos_mir_viejo_round <- matriz_mir(enfermos_round, (as.vector(quantile2mi(0.99741, enfermos_mir_q_round))))
enfermos_gen_viejo_round <- matriz_genes(enfermos_round, (as.vector(quantile2mi(0.99987, enfermos_genes_q_round))))

enfermos_mir_viejo_round[genes,genes] <- enfermos_gen_viejo_round

sanos_mir_viejo_round <- matriz_mir(sanos_round, (as.vector(quantile2mi(0.99741, sanos_mir_q_round))))
sanos_gen_viejo_round <- matriz_genes(sanos_round, (as.vector(quantile2mi(0.99987, sanos_genes_q_round))))

sanos_mir_viejo_round[genes,genes] <- sanos_gen_viejo_round

# viejo round
enfermos_mir_nuevo_round <- matriz_mir(nuevo_enfermos_round, (as.vector(quantile2mi(0.99741, enfermos_mir_q_nuevo_round))))
enfermos_gen_nuevo_round <- matriz_genes(nuevo_enfermos_round, (as.vector(quantile2mi(0.99987, enfermos_genes_q_nuevo_round))))

enfermos_mir_nuevo_round[genes,genes] <- enfermos_gen_nuevo_round

sanos_mir_nuevo_round <- matriz_mir(nuevo_sanos_round, (as.vector(quantile2mi(0.99741, sanos_mir_q_nuevo_round))))
sanos_gen_nuevo_round <- matriz_genes(nuevo_sanos_round, (as.vector(quantile2mi(0.99987, sanos_genes_q_nuevo_round))))

sanos_mir_nuevo_round[genes,genes] <- sanos_gen_nuevo_round



library(minet)
sanos_dpi <- aracne(mim = sanos_mir_viejo, eps = 0.1)
enfermos_dpi <- aracne(mim = enfermos_mir_viejo, eps = 0.1)

sanos_nuevo_dpi <- aracne(mim = sanos_mir_nuevo, eps = 0.1)
enfermos_nuevo_dpi <- aracne(mim = enfermos_mir_nuevo, eps = 0.1)

sanos_round_dpi <- aracne(mim = sanos_mir_viejo_round, eps = 0.1)
enfermos_round_dpi <- aracne(mim = enfermos_mir_viejo_round, eps = 0.1)

sanos_nuevo_round_dpi <- aracne(mim = sanos_mir_nuevo_round, eps = 0.1)
enfermos_nuevo_round_dpi <- aracne(mim = enfermos_mir_nuevo_round, eps = 0.1)


redes <- function(g){
	library(igraph)
	g1 <- graph.adjacency(adjmatrix= g, mode='undirected', diag=F, weighted=T)
	g2 <- get.data.frame(g1)
	return(g2)
}

red_sanos <- redes(sanos_dpi)
red_enfermos <- redes(enfermos_dpi)

red_nuevo_sanos <- redes(sanos_nuevo_dpi)
red_nuevo_enfermos <- redes(enfermos_nuevo_dpi)

red_round_sanos <- redes(sanos_round_dpi)
red_round_enfermos <- redes(enfermos_round_dpi)

red_nuevo_round_sanos <- redes(sanos_nuevo_round_dpi)
red_nuevo_round_enfermos <- redes(enfermos_nuevo_round_dpi)





dim(red_sanos)
dim(red_enfermos)

dim(red_nuevo_sanos)
dim(red_nuevo_enfermos)

dim(red_round_sanos)
dim(red_round_enfermos)

dim(red_nuevo_round_sanos)
dim(red_nuevo_round_enfermos)


> dim(red_sanos)
[1] 33879     3
> dim(red_enfermos)
[1] 29186     3
> 
> dim(red_nuevo_sanos)
[1] 33879     3
> dim(red_nuevo_enfermos)
[1] 29186     3
> 
> dim(red_round_sanos)
[1] 33879     3
> dim(red_round_enfermos)
[1] 29186     3
> 
> dim(red_nuevo_round_sanos)
[1] 33879     3
> dim(red_nuevo_round_enfermos)
[1] 29186     3

> identical(round(red_sanos[,3], digits=6),round(red_nuevo_sanos[,3], digits=6))
[1] TRUE

#------------------------------------------------


enfermos_mir_q_sn <- mi_mirnas(enfermos_sn)
enfermos_genes_q_sn <- mi_genes(enfermos_sn)

sanos_mir_q_sn <- mi_mirnas(sanos_sn)
sanos_genes_q_sn <- mi_genes(sanos_sn)

# matriz sn anterior
enfermos_mir_viejo_sn <- matriz_mir(enfermos_sn, (as.vector(quantile2mi(0.99741, enfermos_mir_q_sn))))
enfermos_gen_viejo_sn <- matriz_genes(enfermos_sn, (as.vector(quantile2mi(0.99987, enfermos_genes_q_sn))))

enfermos_mir_viejo_sn[genes,genes] <- enfermos_gen_viejo_sn

sanos_mir_viejo_sn <- matriz_mir(sanos_sn, (as.vector(quantile2mi(0.99741, sanos_mir_q_sn))))
sanos_gen_viejo_sn <- matriz_genes(sanos_sn, (as.vector(quantile2mi(0.99987, sanos_genes_q_sn))))

sanos_mir_viejo_sn[genes,genes] <- sanos_gen_viejo_sn


sanos_dpi_sn <- aracne(mim = sanos_mir_viejo_sn, eps = 0.1)
enfermos_dpi_sn <- aracne(mim = enfermos_mir_viejo_sn, eps = 0.1)


redes <- function(g){
	library(igraph)
	g1 <- graph.adjacency(adjmatrix= g, mode='undirected', diag=F, weighted=T)
	g2 <- get.data.frame(g1)
	return(g2)
}

red_sanos_sn <- redes(sanos_dpi_sn)
red_enfermos_sn <- redes(enfermos_dpi_sn)


dim(red_sanos)
dim(red_enfermos)

dim(red_sanos_sn)
dim(red_enfermos_sn)

> dim(red_sanos)
[1] 33879     3
> dim(red_enfermos)
[1] 29186     3
> 
> dim(red_sanos_sn)
[1] 38701     3
> dim(red_enfermos_sn)
[1] 36521     3



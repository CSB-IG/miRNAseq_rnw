# load matrix

matrix <- data.matrix(read.table("file.txt"))
colnames(matrix) <- rownames(matrix)
dim(matrix)
isSymmetric(matrix)

# get genes and miRNAs from matrix

mir <- rownames(matrix)[grepl("hsa*", rownames(matrix))]
genes <- rownames(matrix)[grepl("^(?!hsa*)", rownames(matrix), perl=T)]

num_int_genes <- ((length(genes)**2)-length(genes))/2
num_int_mir <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 

# choose a threshold for genes and miRNA edges

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


matrix_mir_q <- mi_mirnas(matrix)
matrix_genes_q <- mi_genes(matrix)

# mismo numero int en genes que int miRNAs
matrix_mir <- matriz_mir(matrix, (as.vector(quantile2mi(0.99741, matrixs_mir_q))))
matrixr_gen <- matriz_genes(matrix, (as.vector(quantile2mi(0.99978, matrix_genes_q))))

matrix_mir[genes,genes] <- matrix_gen


# guardar

write.table(matrix_mir, file="matrix.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

# dpi
library(minet)

matrix_minet <- aracne(mim = matrix_mir, eps = 0.1)

# save

write.table(matrix_minet, file="matrix_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

# save as network

redes <- function(g){
	library(igraph)
	g1 <- graph.adjacency(adjmatrix= g, mode='undirected', diag=F, weighted=T)
	g2 <- get.data.frame(g1)
	return(g2)
}

red_matrix <- redes(matrix_minet)

write.table(red_matrix, file="red_matrix_dpi.txt", sep="\t", col.names=NA, row.names=FALSE, quote=FALSE)




# dpi
networks <- function(prunned_mi_matrix,tol){
	library(minet)
	library(igraph)
	dpi_matrix = aracne(mim = prunned_mi_matrix, eps = tol)
	graph = graph.adjacency(adjmatrix= dpi_matrix, mode='undirected', diag=F, weighted=T)
	edgelist <- get.data.frame(graph)
	return(list(dpi_matrix=dpi_matrix, graph=graph, edgelist=edgelist))
}

network <- networks(prunned_mi_matrix, tol)

# save
write.table(network$dpi_matrix("Names"=rownames(network$dpi_matrix),network$dpi_matrix), file="network_matrix_dpi.txt", row.names=FALSE, sep="\t", quote=F)
write.table(network$edgelist, file="network_edgelist.txt", sep="\t", header=TRUE, row.names=FALSE, quote=FALSE)


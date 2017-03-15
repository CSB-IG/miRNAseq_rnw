
# read p-value = 1 MI squared matrix
mi_matrix <- read.table(file = "MI_matrix.txt", header = T, row.names=1, sep = '\t')
mi_matrix <- data.matrix(frame = mi_matrix, rownames.force = T)
colnames(mi_matrix) <- rownames(mi_matrix)

mir <- rownames(enfermos)[grepl("hsa*", rownames(enfermos))]
genes <- rownames(enfermos)[grepl("^(?!hsa*)", rownames(enfermos), perl=T)]

genes_edges <- ((length(genes)**2)-length(genes))/2
mir_edges <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 

# choose quantile threshold for miRNA and gene edges 
# you may choose an fixed edge number and obtain the corresponding quantile
# number_of_edges = the desired number of edges in the network
# (100/genes|mir_edges)*number_of_edges (before DPI)


matrix_complete <- function(mi_matrix,quantile_mir,quantile_genes,mir,genes){
	# obtain vectors to calculate quantiles
	vector_mir <- c((mi_matrix[mir,mir][upper.tri(mi_matrix[mir, mir])]), as.vector(mi_matrix[mir, genes]))
	vector_genes <- (mi_matrix[genes,genes][upper.tri(mi_matrix[genes, genes])])))
	# convert quantile to mi value
	threshold_mir = as.vector(quantile(vector_mir, probs=quantile_mir, na.rm=T))
	threshold_gen = as.vector(quantile(vector_genes, probs=quantile_genes, na.rm=T))
	# create matrices
	gene_matrix <- mi_matrix
	gene_matrix = gene_matrix[genes,genes]
	gene_matrix[gene_matrix < threshold_gen] <- 0
	mir_matrix <- mi_matrix
	mir_matrix[genes,genes] <- 0
	mir_matrix[mir_matrix < threshold_mir] <- 0
	mir_matrix[genes,genes] <- gene_matrix
	matrix_pval <- mir_matrix
	return(matrix_pval)
}

# quantile_mir = threshold quantile for miR, quantile_genes = threshold quantile genes
prunned_matrix <- matrix_complete(mi_matrix,quantile_mir,quantile_genes,mir,genes)

# save mir gene prunned matrix
write.table(prunned_matrix("Names"=rownames(prunned_matrix),prunned_matrix), file="matrix_mir_gene_pval_threshold.txt", sep="\t", row.names=FALSE, quote=FALSE)



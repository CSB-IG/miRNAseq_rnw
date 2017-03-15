# read squared matrix
tumour <- data.matrix(read.table("matrix_MI_tumour.txt", header=T, row.names=1, sep = '\t'),rownames.force = T)
control <- data.matrix(read.table("matrix_MI_control.txt", header=T, row.names=1, sep = '\t'),rownames.force = T)

colnames(tumour) <- rownames(tumour)
colnames(control) <- rownames(control)

qqplot <- function(tumour, control){
	vector_tumour <- log2(as.vector(tumour[upper.tri(tumour)])+1)
	vector_control <- log2(as.vector(control[upper.tri(control)])+1)
	return(list(vector_tumour=vector_tumour, vector_control=vector_control))
}

qqplot = qqplot(tumour,control)
png(filename = "qqplot.png",
    width = 600, height = 600, units = "px",
     bg = "white",  res = 100)
qqplot(qqplot$vector_tumour, qqplot$vector_control)
dev.off()

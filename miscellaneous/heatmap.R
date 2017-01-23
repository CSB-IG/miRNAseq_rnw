

x <- read.table("matriz_completa.txt", row.names=1, header=T)
x <- data.matrix(x)
mirna <- rownames(x)[grep("^hsa-mir", rownames(x))]
x_mir <- x[mirna,]

library(gplots)
x_mir <- log2(x_mir+1)
heatmap.2(x_mir, col=greenred)

mir200 <- rownames(x)[grep("^hsa-mir-200", rownames(x))]                                                                                                                    
> mir200
[1] "hsa-mir-200a.MIMAT0000682" "hsa-mir-200a.MIMAT0001620"
[3] "hsa-mir-200b.MIMAT0000318" "hsa-mir-200b.MIMAT0004571"
[5] "hsa-mir-200c.MIMAT0000617" "hsa-mir-200c.MIMAT0004657"
> mir200 <- rownames(x)[grep("^hsa-mir-141", rownames(x))]                                                                                                                    
> mir200 <- rownames(x)[grep("^hsa-mir-200", rownames(x))]
> mir141 <- rownames(x)[grep("^hsa-mir-141", rownames(x))]                                                                                                                    
> mir141
[1] "hsa-mir-141.MIMAT0000432" "hsa-mir-141.MIMAT0004598"
> mir429 <- rownames(x)[grep("^hsa-mir-429", rownames(x))]
> mir429
[1] "hsa-mir-429.MIMAT0001536"
> mirna_mir200 <- c(mir200, mir141, mir429)
> mirna_mir200
[1] "hsa-mir-200a.MIMAT0000682" "hsa-mir-200a.MIMAT0001620"
[3] "hsa-mir-200b.MIMAT0000318" "hsa-mir-200b.MIMAT0004571"
[5] "hsa-mir-200c.MIMAT0000617" "hsa-mir-200c.MIMAT0004657"
[7] "hsa-mir-141.MIMAT0000432"  "hsa-mir-141.MIMAT0004598" 
[9] "hsa-mir-429.MIMAT0001536" 


x_mir_200 <- x[mirna_mir200,]
x_mir_200 <- log2(x_mir_200+1)
heatmap.2(x_mir_200, col=greenred)

clusters <- kmeans(t(x_mir_200), 2)

lista <- as.vector(names(sort(clusters$cluster)))

x_mir_200 <- x_mir_200[, lista]
heatmap.2(x_mir_200, col=greenred, )

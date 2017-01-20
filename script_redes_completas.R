# OBSOLETE see Prueba and subtypes_p1_networks
# sin ordenar x
# solo agregar nombres de columnas

> x[is.na(x)] <- 0
> class(x) <- "numeric"

a <- x[15743:16375, 15743:16375]
> dim(a)
[1] 633 633
> sum(is.na(a))
[1] 400689
> 633**2
[1] 400689

# b no tiene miRNAS en rownames

b <- x[1:15742, 1:15742]

adjs_s <- list.files(path="~/TCGA_miRNA_BC/miRNAs_resultados/miRNAS_86_resultados/RNAseq_controles_resultados_ARACNE", pattern="*_1.adj")
> length(adjs_s)
[1] 16225
> length(intersect(adjs, adjs_s))
[1] 15769
genes <- (intersect(adjs, adjs_s))
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
# m son genes que pertenecen a genes pero que no se encuentran en dif
# > 15136-625
# [1] 14511
# > 14511+633
# [1] 15144
# > length(setdiff(faltan, genes))                                                                                                      
# [1] 0
# > length(setdiff(genes, faltan))                                                                                                      
# [1] 14511
# > genes_faltan <- (setdiff(genes, faltan))                                                                                            
# > final <- c(genes_faltan, mir)                                                                                                       
# > length(final)
# [1] 15144
# > 15144+625
# [1] 15769

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


a2 <- x[15743:16375,]

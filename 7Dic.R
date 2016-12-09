


num_int_mir <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 
num_int_genes <- ((length(genes)**2)-length(genes))/2

p= 0.001                                                                                                          │··········································································································································································································
p2mi((p/num_int_mir), num_pacientes_basal)

mi_mirnas <- function(x){
	return(c((x[mir,mir][upper.tri(x[mir, mir])]), as.vector(x[mir, genes])))
}

mi_genes <- function(x){
	return(x[genes,genes][upper.tri(x[genes, genes])])
}


basal_mir_q <- mi_mirnas(basal)
her2_mir_q <- mi_mirnas(her2)
luma_mir_q <- mi_mirnas(luma)
lumb_mir_q <- mi_mirnas(lumb)
casos_mir_q <- mi_mirnas(casos)
sanos_mir_q <- mi_mirnas(sanos)

basal_genes_q <- mi_genes(basal)
her2_genes_q <- mi_genes(her2)
luma_genes_q <- mi_genes(luma)
lumb_genes_q <- mi_genes(lumb)
casos_genes_q <- mi_genes(casos)
sanos_genes_q <- mi_genes(sanos)

# quantile(her2_genes_q, probs=seq(0.995, 1, 0.00001), na.rm=T)

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


basal_mir_p <- matriz_mir(basal, (as.vector(quantile2mi(0.99770, basal_mir_q))))
her2_mir_p <- matriz_mir(her2, (as.vector(quantile2mi(0.99770, her2_mir_q))))
luma_mir_p <- matriz_mir(luma, (as.vector(quantile2mi(0.99770, luma_mir_q))))
lumb_mir_p <- matriz_mir(lumb, (as.vector(quantile2mi(0.99770, lumb_mir_q))))
casos_mir_p <- matriz_mir(casos, (as.vector(quantile2mi(0.99770, casos_mir_q))))
sanos_mir_p <- matriz_mir(sanos, (as.vector(quantile2mi(0.99770, sanos_mir_q))))

basal_gen_p <- matriz_genes(basal, (as.vector(quantile2mi(0.99987, basal_genes_q))))
her2_gen_p <- matriz_genes(her2, (as.vector(quantile2mi(0.99987, her2_genes_q))))
luma_gen_p <- matriz_genes(luma, (as.vector(quantile2mi(0.99987, luma_genes_q))))
lumb_gen_p <- matriz_genes(lumb, (as.vector(quantile2mi(0.99987, lumb_genes_q))))
casos_gen_p <- matriz_genes(casos, (as.vector(quantile2mi(0.99987, casos_genes_q))))
sanos_gen_p <- matriz_genes(sanos, (as.vector(quantile2mi(0.99987, sanos_genes_q))))

basal_p <- basal_mir_p
her2_p <- her2_mir_p
luma_p <- luma_mir_p
lumb_p <- lumb_mir_p
casos_p <- casos_mir_p
sanos_p <- sanos_mir_p

basal_p[genes,genes] <- basal_gen_p
her2_p[genes,genes] <- her2_gen_p
luma_p[genes,genes] <- luma_gen_p
lumb_p[genes,genes] <- lumb_gen_p
casos_p[genes,genes] <- casos_gen_p
sanos_p[genes,genes] <- sanos_gen_p

write.table(basal_p, file="basal_prunning.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(her2_p, file="her2_prunning.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(luma_p, file="luma_prunning.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(lumb_p, file="lumb_prunning.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(casos_p, file="casos_prunning.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_p, file="sanos_prunning.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

basal_minet <- aracne(mim = basal_p, eps = 0.1)
her2_minet <- aracne(mim = her2_p, eps = 0.1)
luma_minet <- aracne(mim = luma_p, eps = 0.1)
lumb_minet <- aracne(mim = lumb_p, eps = 0.1)
casos_minet <- aracne(mim = casos_p, eps = 0.1)
sanos_minet <- aracne(mim = sanos_p, eps = 0.1)

write.table(basal_minet, file="basal_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(her2_minet, file="her2_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(luma_minet, file="luma_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(lumb_minet, file="lumb_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(casos_minet, file="casos_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(sanos_minet, file="sanos_dpi.txt", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

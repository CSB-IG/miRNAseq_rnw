
# 752
genes <- read.table("RNAseq_752_norm.rsem.aracne.txt", sep="\t", header=T, row.names=1)

miRNAs <- read.table("miRNA_752_normTMM.txt", sep="\t", header=T, row.names=1)

identical(colnames(genes), colnames(miRNAs))
[1] TRUE

> red <- rbind(miRNAs, genes)
dim(red)

write.table(red, quote = F, sep = "\t", row.names = T, col.names = NA, file = "matriz_completa.txt")

### 86
genes_enfermos <- read.table("enfermos.86.txt", sep="\t", header=T, row.names=1)
genes_sanos <- read.table("sanos.86.txt", sep="\t", header=T, row.names=1)
enfermos <- read.table("miRNAs_86_maduros_norm_enfermos.txt", sep="\t", header=T, row.names=1)
sanos <- read.table("miRNAs_86_maduros_norm_sanos.txt", sep="\t", header=T, row.names=1)

identical(colnames(genes_enfermos), colnames(enfermos))
identical(colnames(genes_sanos), colnames(sanos))

length(intersect(rownames(genes_enfermos), rownames(genes_sanos)))
[1] 15136

genes_enfermos <- genes_enfermos[intersect(rownames(genes_enfermos), rownames(genes_sanos)), ]
genes_sanos <- genes_sanos[intersect(rownames(genes_enfermos), rownames(genes_sanos)), ]
identical(rownames(genes_enfermos), rownames(genes_sanos))
identical(gsub('.{3}$', '',colnames(genes_enfermos)), gsub('.{3}$', '',colnames(genes_sanos)))

red_enfermos <- rbind(enfermos, genes_enfermos)
red_sanos <- rbind(sanos, genes_sanos)
dim(red_enfermos)
dim(red_sanos)

write.table(red_enfermos, quote = F, sep = "\t", row.names = T, col.names = NA, file = "matriz_completa_enfermos.txt")
write.table(red_sanos, quote = F, sep = "\t", row.names = T, col.names = NA, file = "matriz_completa_sanos.txt")




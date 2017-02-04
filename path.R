
enfermos <- read.table("RNAseq_752_norm.rsem.aracne.txt", )
dim(enfermos)

sanos <- read.table("sanos.86.txt")
dim(sanos)

genes <- intersect(rownames(enfermos), rownames(sanos))

enfermos <- enfermos[genes,]
sanos <- sanos[genes,]

enf <- rep("0", 752)
san <- rep("1", 86)
enfermos_m <- rbind(enf,enfermos)
sanos_m <- rbind(san, sanos)

matriz <- cbind(enfermos_m, sanos_m)

id <- read.table("annot.txt")
id <- id[rownames(enfermos),]

nombres <- c("NORMALS", id[,2])
geneid <- c("NORMALS", id[,1])

rownames(matriz) <- nombres
write.table(matriz, "matriz_path_ensbl.txt", sep="\t", quote=F)

rownames(matriz) <- geneid
write.table(matriz, "matriz_path_geneid.txt", sep="\t", quote=F)


# add names in the first colname

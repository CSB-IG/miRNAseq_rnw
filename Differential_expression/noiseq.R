
library(plyr)
library(NOISeq)


# read expression matrices that were used for mi computation
enfermos <- read.table(file = "enfermos.86.txt", header = T, sep = '\t')
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')

# create a vector with gene IDs from both matrices
enfermos <- as.vector(enfermos[ ,1])
sanos <- as.vector(sanos[ ,1])


# gene list from genes in common
genes <- sort(intersect(enfermos, sanos))


# read expression matrices before filtering 
enfermos <- read.csv("RNAseq_norm_casos_86.csv")
sanos <- read.csv("RNAseq_norm_control_86.csv")

# only use genes with GENE ID
enfermos <- enfermos[30:20531, ]
sanos <- sanos[30:20531, ]
dim(enfermos)
dim(sanos)

# assign gene IDs
ID <- scan(file="GeneID", skip=29, what="")
row.names(enfermos) <- ID
row.names(sanos) <- ID

# filter matrix by a list: matrix[list, ]
enfermos <- enfermos[genes, ]
sanos <- sanos[genes, ]

# create a single matrix with expression data from cancer and control tissue
RNAseq <- cbind(enfermos, sanos)
write.table(RNAseq, quote = F, sep = "\t", row.names = T, col.names = NA, file = "noiseq.txt")

# factors in the same order that cbind: enfermos -> sanos
myfactors <- data.frame(Pacientes = c( rep('enfermos', 86), rep('sanos', 86)))

# creating NOIseq object
x <- readData(data= RNAseq, factors=myfactors)

str(x)
head(assayData(x)$exprs)
head(pData(x))

# QC report
QCreport(x, samples = NULL, factor = "Pacientes", norm = FALSE)

# differential expression
exp <- noiseq(x, factor="Pacientes", k=NULL, norm="n", pnr= 0.2, nss= 5, v=0.02, lc=1, replicates = "no") 
head(exp@results[[1]])

# differential expressed genes
exp.deg = degenes(exp, q = 0.7, M = NULL)

# up-regulated in first condition
exp.deg1 = degenes(exp, q = 0.7, M = "up")

# down-regulated in first condition
exp.deg1 = degenes(exp, q = 0.7, M = "down")

# plot
pdf("exp_plot.pdf",width=14,height=10)
DE.plot(exp, q=0.7, graphic="expr", log.scale=TRUE)
dev.off()
pdf("MD_plot.pdf",width=14,height=10)
DE.plot(exp, q=0.7, graphic="MD")
dev.off()

# save results
write.table(exp.deg, file = "dif_exp_RNAseq.txt", sep = "\t", col.names= T, row.names= T, quote = F )





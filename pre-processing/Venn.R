library(gplots)

datos<-read.csv('suma_miRNAhiseq-only.csv')
RNAseq_casos1 <- datos[,1]
RNAseq_controles1 <- datos[,3]
miRNAseq_casos1 <- datos[,2]
miRNAseq_controles1 <- datos[,4]

#remove 4 last digits from the TCGA ID so they are comparable
#between controls and tumors
RNAseq_casos <- gsub('.{4}$', '', RNAseq_casos1)
RNAseq_controles <- gsub('.{4}$', '', RNAseq_controles1)
miRNAseq_casos <- gsub('.{4}$', '', miRNAseq_casos1)
miRNAseq_controles <- gsub('.{4}$', '', miRNAseq_controles1)

#removes empty values
x1 <- (RNAseq_casos)
x2 <- head(RNAseq_controles,-982)
x3 <- head(miRNAseq_casos,-340)
x4 <- head(miRNAseq_controles,-1008)

#venn diagram
venn(list(mRNA_casos=x1, mRNA_controles=x2, miRNA_casos=x3, miRNA_controles=x4))



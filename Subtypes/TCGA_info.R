
abrir txt TCGA en R
TCGA <- read.table(file="TCGA_BRCA.txt", sep="\t", header=T, stringsAsFactors=F)
TCGA <- TCGA[,1:5]
TCGA[TCGA$Tissue == 'Metastatic',]
TCGA <- TCGA[-(364:366),]
pacientes <- TCGA$bcr_patient_barcode
Pacientes <- sub('.*(?=.{4}$)','', pacientes, perl=T)
setwd("/home/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/")
enfermos752 <- read.table(file = "RNAseq_752_norm.rsem.aracne.txt", header = T, sep = '\t')
muestras <- colnames(enfermos752[-1])
Muestras <- gsub('.{3}$', '', muestras)
library(gplots)
pdf("subtipo.pdf",width=14,height=10)
venn(list(Pacientes, Muestras))
dev.off()
Sub <- intersect(Pacientes, Muestras)
rownames(TCGA) <- Pacientes
Subtipo <- TCGA[Sub,]
write.table(Subtipo, quote = F, sep = "\t", row.names = T, col.names = NA, file = "Subtipo_170_de_752.txt")
sum(Subtipo$PAM50 == "Basal")
33
sum(Subtipo$PAM50 == "LumA")
81
sum(Subtipo$PAM50 == "LumB")
33
sum(Subtipo$PAM50 == "Her2")
22

### VENN PLOT
TCGA <- read.table("Subtipos_170_RNAseq.txt", header=T, sep="\t")
genefu <- read.table("Subtipos_752_RNAseq.txt", header=T, sep="\t")

LAg <- genefu[genefu$PAM50 == "LumA",1]
LBg <- genefu[genefu$PAM50 == "LumB",1]
H2g <- genefu[genefu$PAM50 == "Her2",1]
Bg <- genefu[genefu$PAM50 == "Basal",1]
Ng <- genefu[genefu$PAM50 == "Normal",1]

LAt <- TCGA[TCGA$PAM50 == "LumA",1]
LBt <- TCGA[TCGA$PAM50 == "LumB",1]
H2t <- TCGA[TCGA$PAM50 == "Her2",1]
Bt <-  TCGA[TCGA$PAM50 == "Basal",1]
Nt <- TCGA[TCGA$PAM50 == "Normal",1]

pdf("PAM50vs.pdf",width=10,height=10)
venn(list(TCGA=TCGA[,1], genefu=genefu[,1]))
venn(list(TCGA_LumA=LAt, genefu=LAg))
venn(list(TCGA_LumB=LBt, genefu=LBg))
venn(list(TCGA_Her2=H2t, genefu=H2g))
venn(list(TCGA_Basal=Bt, genefu=Bg))
venn(list(TCGA_Normal=Nt, genefu=Ng))
dev.off()

genefu[genefu$Paciente == "A13Z.01",]

her2 <- c("A18P.01", "A12Q.01", "A0CL.01", "A138.01", "A0RH.01", "A14V.01", "A0B7.01", "A12T.01", "A130.01", "A1EV.01")
6 LumA/ 4 LumB














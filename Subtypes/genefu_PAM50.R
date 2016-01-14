
awk -F'|' '{ print $1"\t"$1"\t"$2 }' Gene_ID | less -S

library(genefu)
# genefu input: normalized expression matrix with geneID as columns & samples as rows
# 		annotation table: at least one col = probe (GeneID) & EntrezGene.ID (EntezID)

dataset <- read.table("RNAseq_752_norm.rsem.aracne.txt", header=T, sep="\t", row.names=1)
dataset <- as.data.frame(t(dataset))
dannot <- read.table("annot.txt", header=F, sep="\t", row.names=1)
colnames(dannot)<-c("probe", "EntrezGene.ID")


# pam50.scale
# pam50.robust
PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50, data=dataset, annot=dannot, do.mapping=TRUE, verbose=TRUE)
table(PAM50Preds$subtype)

LumA <- names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumA")]
LumB <- names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumB")]
Basal <- names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "Basal")]
Her2 <- names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "Her2")]
Normal_like <- names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "Normal")]

write.table(PAM50Preds$subtype, quote = F, sep = "\t", row.names = T, file = "Subtipos_pam50.scale_752_RNAseq.txt")

# https://www.biostars.org/p/102212/
# First, make sure you have an annotation dataset. They can be found in bioconductor (i.e. annot.nkis, org.Hs.egALIAS2EG). 
# Then make sure that, if not already, this annotation dataset has the column name "EntrezGene.ID" 
# plus whatever other columns in the dataset present such as gene_name. 
# Now, if you are doing predictions (such as for PAM50), make sure you are using intrinsic.cluster.predict() 
# and not intrinsic.cluster(). At the end you should have your 3 objects for your predictions: 
# the pam50 model found in genefu (data(pam50)), your matrix that SHOULD HAVE SAMPLES IN ROWS AND GENES IN COLUMNS 
# (this seems to be the most confusing part) and lastly the annotation data.frame with a column name "EntrezGene.ID."

# > table(PAM50Preds$subtype)

# Basal   Her2   LumA   LumB Normal 
#   150    109    305    124     64 





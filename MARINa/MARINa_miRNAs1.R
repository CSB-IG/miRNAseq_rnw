######## LOADS ########

setwd("~/miRNAseq_rnw/MARINa/")
suppressPackageStartupMessages(library(mixtools))
library(ssmarina)
load("exp_set.RData")


# Regulon generation
adjfile <- file.path("enfermos_reguls_network.adj")
regulon <- aracne2regulon(adjfile, eset)

save(regulon, file = "./regulon_miRNAs1.Rdata")

########  S I G N A T U R E ######## 

# Row by row T test comparation of the two matrix
signature <- rowTtest(eset[, case], eset[, control])

# z-score values for the GES
signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]


########  N U L L  M O D E L ######## 

# NULL model by sample permutations
nullmodel <- ttestNull(eset[, case], eset[, control], per=1000, repos=T)


save(signature,nullmodel,file="./sign_null_miRNAs1.Rdata")


to_clean <- ls()

########  M A R I N A ######## 

mrs <- marina(signature, regulon, nullmodel)

mrs_noledges <- mrs

save(mrs,top100mrs,mrs_noledges,file="miRNAs1_MARINa.RData")

# Leading-edge analysis
mrs <- ledge(mrs)

# Ploting masters regulators
pdf("Top10mrs_miRNAs1.pdf",width=6,height=7)
plot(mrs_noledges, cex=.7)
dev.off()

print(summary(mrs_noledges))

top100mrs <-summary(mrs,mrs=100)
write.table(top100mrs,file = "./miRNAs1_top100mrs.txt",quote = FALSE, sep = "\t")


#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
mrsshadow <- shadow(mrs, pval=25)

write.table(summary(mrsshadow)$Shadow.pairs, file="miRNAs1_shadow_pairs.sif", quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Synergy analysis
mrs <- marinaCombinatorial(mrs, regulators=25)

mrs <- marinaSynergy(mrs)

# ploting synergy regulators
pdf("Synergy_miRNAs1.pdf",width=11,height=7)
plot(mrs,mrs=10, cex=.7)
dev.off()

######### F I N A L  S A V E  ##########
rm(list=to_clean,to_clean) #cleaning
save.image("miRNAs1_MARINa.RData")

print("All done!")
controles <- read.table("controles_int.txt")
casos <- read.table("casos_int.txt")

hcasos <- hist(sort(casos[,1]), 
	las=0.1, 
	xlim=c(0,1),
	ylim=c(0,400),
	breaks=1181)
hcontroles <- hist(sort(controles[,1]), 
	xlab="Mutual Information", 
	ylab="Frequency",
	las=0.1, 
	xlim=c(0,1),
	ylim=c(0,200),
	breaks=611)
	
hcasos$counts[hcasos$counts==0] <- NA
hcontroles$counts[hcontroles$counts==0] <- NA

pdf("prueba.pdf", width=15, height=10)
plot(hcasos$mids, hcasos$counts,
	pch=20,
	xlab= NA,
	ylab=NA,
	xlim=c(0,1),
	ylim=c(0,200),
	col="red")
points(hcontroles$mids, hcontroles$counts,
	pch=20,
	col="black")
dev.off()



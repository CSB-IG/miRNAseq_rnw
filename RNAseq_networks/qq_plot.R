# read squared matrix
we <- read.table(file = "enfermos_adj.txt", header = T, sep = '\t')
ws <- read.table(file = "sanos_adj.txt", header = T, sep = '\t')

colnames(we) <- rownames(we)
colnames(ws) <- rownames(ws)

we <- data.matrix(frame = we, rownames.force = T)
ws <- data.matrix(frame = ws, rownames.force = T)

we[we == 0] <- NA
ws[ws == 0] <- NA

we <- log10(we)
ws <- log10(ws)

# verify that is symmetric
isSymmetric(we)
isSymmetric(ws)
we[upper.tri(we, diag=T)] <- NA
ws[upper.tri(ws, diag=T)] <- NA

pdf("qq_plot_es.pdf",width=5,height=5)
qqplot(we, ws)
dev.off()

png(filename = "qqplot86.png",
    width = 600, height = 600, units = "px",
     bg = "white",  res = 100)
qqplot(we, ws)
dev.off()
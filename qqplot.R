we <- read.table(file = "enfermos_p1_full_norm_sif.txt", header = T, sep = '\t')
ws <- read.table(file = "sanos_p1_full_norm_sif.txt", header = T, sep = '\t')

we <- as.numeric(we[,3])
ws <- as.numeric(ws[,3])

we[we == 0] <- NA
ws[ws == 0] <- NA

png(filename = "qqplot86.png",
    width = 600, height = 600, units = "px",
     bg = "white",  res = 100)
qqplot(we, ws)
dev.off()


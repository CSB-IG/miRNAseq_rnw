# data set
# mutual information values
mi <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# data with or without dpi in order
dpi <- c('sin.dpi', 'sin.dpi', 'sin.dpi', 'sin.dpi','sin.dpi', 'sin.dpi', 'sin.dpi','sin.dpi','sin.dpi', 'con.dpi', 'con.dpi', 'con.dpi', 'con.dpi','con.dpi', 'con.dpi', 'con.dpi', 'con.dpi', 'con.dpi')

# log 10 of the number of interaction for that mi value interaction
log10_interacciones <- log10(c(127465765, 8380507, 1160687, 211581, 51689, 15320, 4560, 1198, 274, 126871142, 8159284, 1111557, 200424, 49708, 15147, 4557, 1198, 274))

# create data.frame
w <- cbind(mi,log10_interacciones, dpi)

# graph
library(ggplot2)
w <- as.data.frame(w)
w$log10_interacciones <- as.numeric(as.character(w$log10_interacciones))
ggplot(w, aes(mi, log10_interacciones, fill = dpi)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = c("maroon3", "deeppink1"))

# line plot
# mutual information values (without dpi)
mi_sindpi <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# log 10 of the number of interactions (without dpi) 
log10_interacciones_sindpi<- log10(c(127465765, 8380507, 1160687, 211581, 51689, 15320, 4560, 1198, 274))

# get data.frame (without dpi)
w_sindpi <- cbind(mi_sindpi, log10_interacciones_sindpi)
plot(w_sindpi, xlab = "MI" ,ylab = "log10_interacciones", type = "o", pch= 21, bg="#dd1c77")
par(new=TRUE)

#  mutual information values (with dpi)
mi_dpi <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# log 10 of the number of interactions (with dpi)
log10_interacciones_dpi <- log10(c(126871142, 8159284, 1111557, 200424, 49708, 15147, 4557, 1198, 274))

# get data.frame (with dpi)
w_dpi <- cbind(mi_dpi,log10_interacciones_dpi)

# plot and facet by categorie
plot(w_dpi, xlab = "MI" ,ylab = "log10_interacciones", type = "o", pch= 21, bg="pink")

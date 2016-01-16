# mutual information (mi) values
MI <-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8) 

# corresponding p-value for mi
p_value <- -log10(c(4.369537e-23, 6.601615e-46, 9.973897e-69, 1.506883e-91, 2.27664e-114, 3.439609e-137, 5.196654e-160, 7.851244e-183))

# plot
plot(MI, p_value)

# plot with colors
plot(x = MI, y = p_value, xlab = "MI" ,ylab = "-log10_pvalue", type = "o", col.lab="black", pch= 21, bg="#dd1c77", col="#980043")


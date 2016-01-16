# p-value to mutual information converter 

p2mi <- function(x, n){
  alfa <- 1.062
  beta <- -48.7
  gamma <- -0.634
  mi <- (alfa - log(x)) / ((-beta) + (-gamma)*n)
  return(mi)
  }

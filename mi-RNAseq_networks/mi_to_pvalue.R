# mutual information to p-value converter

mi2p <- function(mi, n){
  alfa <- 1.062
  beta <- -48.7
  gamma <- -0.634
  x <- exp(-((mi*((-beta) + (-gamma)*n)) - alfa))
  return(x)
}


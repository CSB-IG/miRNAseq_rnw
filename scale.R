scale_mtx <-function(m){
	if(min(m[upper.tri(m)]) == 0){
		a <- (m-min(m))/(max(m)-min(m))
		return(a)
	}else{
		a <- (m - min(m[upper.tri(m)]))/(max(m)-min(m[upper.tri(m)]))
		diag(a) <- 0
		return(a)
		}
}

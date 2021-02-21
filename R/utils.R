# Elbow knee --------------------------------------------------------------

#' elbow_knee
#'
#'Identifies the principle components (PC) to use by using elbow knee method
#'
#' @param x float: the standard deviation of each principle component
#' @param addpca int: the PCs to add to those considered
#'
#' @return int: Number of PCs to consider
#' @export
elbow_knee = function(x, addpca = 5){
  y = c(x[1],1)
  z = c(x[length(x)],length(x))
  dist.vec = c()
  for(i in 1:length(x)){
    i.pt = c(x[i],i)
    dist.vec[i] = dist.pt(i.pt,y,z)
  }
  return(which.max(dist.vec)+addpca)
}

# Distance between line and point --------------------------------------------------------------

#' dist.pt
#'
#' Distance between a line and point
#'
#' @param x int: the point
#' @param y int: the standard deviation of the 1st PC
#' @param z int: the standard deviation of the last PC
#'
#' @return
#' @export
dist.pt = function(x,y,z){
  v1 <- y - z
  v2 <- x - y
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
}


# Progress bar with text --------------------------------------------------

#' pog.txt
#'
#' @param i index
#' @param n total number of iteration
#'
#' @return
#' @export
prog.txt = function(i,n){
  cat("\r", paste0(round(100*i/n,2), "% complete"))
}

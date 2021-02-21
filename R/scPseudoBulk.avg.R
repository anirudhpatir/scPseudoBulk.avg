
# scPseudoBulk.avg --------------------------------------------------------

#' scPseudoBulk.avg
#'
#' @param norm data, matrix/dataframe: Normalized count matrix (rows = genes, columns = cells)
#' @param clusters data, factors: vector of clusters for each cell
#' @param th.exp.pct threshold, int: quantile cap of expression values within a cluster to remove spikes in expression (default = 90)
#' @param th.drp.pct threshold, int: minimum percent of dropout cells within a cell cluster below which expression values are set to 0 (default = 10)
#' @param th.drp.cell threshold, int: minimum number of dropout cells within a cell cluster below which expression values are set to 0 (default = 3)
#'
#' @return
#' @export
scPseudoBulk.avg = function(
  norm = NULL,
  clusters = NULL,
  th.exp.pct = 90,
  th.drp.pct = 10,
  th.drp.cell = 3
  ){

  #### Filter, Cap and average matrix ####

  norm.flt = list()
  for(i in 1:length(levels(clusters))){
    prog.txt(i,length(levels(clusters)))
    tmp.ind = clusters==levels(clusters)[i]
    tmp.norm = norm[, tmp.ind]
    norm.flt[[i]] = apply(tmp.norm,1,function(x) flt.cluster(x,
                                                       th.exp.pct,
                                                       th.drp.pct,
                                                       th.drp.cell))
  }
  norm.flt = do.call(cbind.data.frame, norm.flt)
  colnames(norm.flt) = levels(clusters)

  #### Remove genes with no expression ####

  kp = apply(norm.flt,1,function(x) any(x>0))
  norm.flt = norm.flt[kp,]

}

# scPseudoBulk.avg --------------------------------------------------------


#' Title
#'
#' @param x data, int: vector of expression values
#' @param th.exp.pct threshold, int: quantile cap of expression values within a cluster to remove spikes in expression (default = 90)
#' @param th.drp.pct threshold, int: minimum percent of dropout cells within a cell cluster below which expression values are set to 0 (default = 10)
#' @param th.drp.cell threshold, int: minimum number of dropout cells within a cell cluster below which expression values are set to 0 (default = 3)
#'
#' @return
#' @export
flt.cluster = function(x, th.exp.pct, th.drp.pct, th.drp.cell){

  #### Filter on dropout #####

  mx.drp = max( round(th.drp.pct*length(x)/100), th.drp.cell)
  if(sum(x>0) <= mx.drp){
    x = rep(0,length(x))
  }

  #### Cap on expression value #####

  if(any(x>0)){
    mx.exp.pct = quantile(x, th.exp.pct/100)
    mx.x = max(x[which(x<mx.exp.pct)])
    if(mx.x == 0){x[x>mx.x] = mx.exp.pct}
    else{x[x>mx.x] = mx.x}
  }

  return(mean(x))

}

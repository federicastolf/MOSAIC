library(vegan)

sourceCpp("helper/Cwrapper.cpp")

################################################################################
#------------------------# ecology dissimilarities #---------------------------#
################################################################################

jaccard_diss = function(x, y) {
  intersection = sum(x & y)
  union = sum(x | y)
  if (union == 0) return(0)  # define dissimilarity as 0 if both vectors are all 0
  return(1 - intersection / union)
}

sorensen_diss = function(x, y) {
  A = sum(x & y)
  sx = sum(x)
  sy = sum(y)
  if ((sx + sy) == 0) return(0)
  sr = (2 * A) / (sx + sy)  # same as 2A/(2A+B+C)
  return(1-sr)
}


Betadiss_data = function(Y, n, method){
  ngroup = length(n)
  # to construct Yq
  ends_n = cumsum(n)
  starts_n = c(1, head(ends_n, -1) + 1)
  # sitelb = levels(lab_group)
  # ngroup = length(sitelb)
  dissAll = rep(0, (ngroup*(ngroup-1))/2)
  id=1
  for(q in 1:(ngroup-1)){
    # gname = sitelb[q]
    # idg = which(lab_group==gname)
    # Yg = Y[idg,]
    Yg = Y[starts_n[q]:ends_n[q], ]
    indq = as.integer(colSums(Yg == 1) > 0)
    for(p in (q+1):ngroup){
      # pname = sitelb[p]
      # idp = which(lab_group==pname)
      # Yp = Y[idp,]
      Yp = Y[starts_n[p]:ends_n[p], ]
      indp = as.integer(colSums(Yp == 1) > 0)
      if(method=="Jaccard"){
        dissAll[id] = jaccard_diss(indp, indq)
        id = id+1
      }
      if(method=="Sorensen"){
        dissAll[id] = sorensen_diss(indp, indq)
        id = id+1
      }
    }
  }
  return(dissAll)
}
APC_correction = function(MI_mat) {
  L = ncol(MI_mat)
  
  APC_mat = matrix(nrow = L, ncol = L)
  meanrow = apply(MI_mat, 1, function(x) {
    sum(x, na.rm = T)
  })
  meanall = sum(MI_mat, na.rm = T)
  
  for (i in 1:L) {
    for (j in i:L) {
      APC_mat[i, j] = meanrow[i] * meanrow[j] / meanall
      APC_mat[j, i] = APC_mat[i, j]
    }
    APC_mat[i, i] = NA
  }
  
  MIp= MI_mat - APC_mat
  return(MIp)
}

APC_correction_tcol=function(mat_tcol,dim){
  mat=matrix(0,dim,dim)
  
  mat[as.matrix(mat_tcol[,c(1,2)])]=mat_tcol[,3]
  mat[as.matrix(mat_tcol[,c(2,1)])]=mat_tcol[,3]
  
  mat_apc=APC_correction(mat)
  return(mat_apc)
}
falign_R=function(score_mtx, rows,cols){
  max_sco = 0;
  # sco[rows+1,cols+1];
  gap_o=-1
  sco=matrix(0, rows+1,cols+1)
  for (i in 2:(rows+1)){
    for (j in 2:(cols+1)){
      A = sco[i-1,j-1] + score_mtx[i-1,j-1]
      D = sco[i-1,j]
      R = sco[i,j-1]

      if (A >= R) {
        if (A >= D) {
          sco[i, j] = A
        } else{
          sco[i,j] = D
        }
      }
      else{
        if (R >= D) {
          sco[i,j] = R
        } else{
          sco[i,j] = D
        }
      }

      if(sco[i,j] > max_sco){max_sco = sco[i,j]}
    }
  }
  return(max_sco)
  # return(sco)
}

align_R=function(score_mtx, gap_open=-1,gap_e=0.2,debug=FALSE){
  rows=nrow(score_mtx)
  cols=ncol(score_mtx)

  sco=matrix(0, rows+1,cols+1)
  label=matrix(0, rows+1,cols+1)
  max_sco=0
  a2b=integer(rows)
  a2b[]=-1
  max_i=0
  max_j=0

  for (i in 2:(rows+1)){

    for (j in 2:(cols+1)){
      A = sco[i-1,j-1] + score_mtx[i-1,j-1]
      D = sco[i-1,j]
      R = sco[i,j-1]
      if(label[i-1,j] == 1){D =D+ gap_open}else{D =D+ gap_open * gap_e}
      if(label[i,j-1] == 1){R =R+ gap_open}else{R =D+ gap_open * gap_e}
      # if(label[i-1,j] == 1){D =D+ gap_b[j-1]}else{D =D+ gap_b[j-1] * gap_e}
      # if(label[i,j-1] == 1){R =R+ gap_a[i-1]}else{R =D+ gap_a[i-1] * gap_e}

      if(A <= 0 && D <= 0 && R <= 0){label[i,j] = 0;sco[i,j] = 0;}
      else{
        if(A >= R){if(A >= D){label[i,j] = 1;sco[i,j] = A;}else{label[i,j] = 2;sco[i,j] = D;}}
        else{if(R >= D){label[i,j] = 3;sco[i,j] = R;}else{label[i,j] = 2;sco[i,j] = D;}}
        if(sco[i,j] > max_sco){max_i = i;max_j = j;max_sco = sco[i,j];}
      }
    }
  }

  i = max_i;
  j = max_j;

  while (1) {
    if (label[i,j] == 0) {
      break
    }
    else if (label[i,j] == 1) {
      a2b[i - 1] = j - 1
      i=i-1
      j=j-1
    }
    else if (label[i,j] == 2) {
      i=i-1
    }
    else if (label[i,j] == 3) {
      j=j-1
    }
  }

  if (debug){
    return(sco)
  }

  return(a2b)

}


#' read_mrf read a MRF model from GREMLIN output.
#'
#' @param filemrf
#'
#' @return list object of mrf
#' @export
#' @import dplyr
#'
read_mrf = function(filemrf) {

  myalphabet = c("a", "u", "c", "g", "-")

  v1 = data.table::fread(cmd = paste("grep '^V'", filemrf))
  names(v1)[-1] = myalphabet

  w1 = data.table::fread(cmd = paste("grep  '^W'", filemrf))


  len=nrow(v1)
  len_a=length(myalphabet)

  w1$i=as.integer(gsub(".*\\[(.*?)\\]\\[(.*?)\\].*","\\1",w1$V1))+1

  w1$j=as.integer(gsub(".*\\[(.*?)\\]\\[(.*?)\\].*","\\2",w1$V1))+1

  array_j = array(0, dim = c(len, len, len_a, len_a))

  for (m in 1:nrow(w1)) {
    id_i = w1$i[m]
    id_j = w1$j[m]

    mat = matrix(as.matrix(w1[m, 2:26]), 5, 5, byrow = TRUE)
    array_j[id_i, id_j, , ] = mat

  }



  ids = expand.grid(1:nrow(v1), 1:nrow(v1)) %>% filter(Var2 > Var1) %>% arrange(Var1)
  w1_mat = as.matrix(w1[, 2:(1 + length(myalphabet) * length(myalphabet))])

  mat_score = sapply(1:nrow(w1_mat), function(i) {
    tmpmat = matrix(w1_mat[i,], 5, 5)
    score = sqrt(sum(tmpmat[-5,-5] ^ 2))
    return(score)
  })
  ids$score = mat_score

  mat_mrf = matrix(0, nrow(v1), nrow(v1))

  mat_mrf[as.matrix(ids[, c(1, 2)])] = ids[, 3]
  mat_mrf[as.matrix(ids[, c(2, 1)])] = ids[, 3]

  mat_apc = APC_correction(mat_mrf)



  mrf = list(
    len = len,
    h = v1,
    j = w1,
    mat_mrf = mat_mrf,
    mat_apc = mat_apc,
    array_j=array_j
  )
  mrf_mat=mrf2mrf_mat(mrf)
  mrfh = as.matrix((mrf$h[, 2:6]))

  mrf$mrf_mat=mrf_mat
  mrf$mrf_h=mrfh

  return(mrf)
}

retrieve_matj=function(i,a,j,b,mat_j,len_a){
  return(mat_j[(i-1)*len_a+a,(j-1)*len_a+b])
}

mrf2mrf_mat = function(mrf) {
  myalphabet = c("a", "u", "c", "g", "-")

  len = mrf$len
  len_a = length(myalphabet)


  mat_j = matrix(0, len * len_a, len * len_a)

  w1 = mrf$j
  for (m in 1:nrow(w1)) {
    id_i = w1$i[m]
    id_j = w1$j[m]

    id_ia = id2_to_id1(1, id_i, len_a)
    id_ja = id2_to_id1(1, id_j, len_a)

    mat = matrix(as.matrix(w1[m, 2:26]), 5, 5, byrow = TRUE)
    # array_j[id_i, id_j, ,] = mat

    mat_j[id_ia:(id_ia + len_a - 1), id_ja:(id_ja + len_a - 1)] = mat
  }
  return(mat_j)
}

id1_to_id2=function(id,dim){
  # j = id %/% (dim +1)+1
  i = id %% (dim)
  j=ceiling(id/dim)
  i[i==0]=dim
  return(data.frame(i=i,j=j))
}

id2_to_id1=function(i,j,dim){
  return((j-1)*dim+i)
}


#' Encode RNA sequence
#'
#' @param seq
#'
#' @return a list
#' @export
#'
encode_seq = function(seq) {
  myalphabet = c("a", "u", "c", "g")
  seq_int = match(seq, table = myalphabet, nomatch = 0) - 1 # 0 based
  seq_int_ungapped = seq_int[seq_int > -1]
  seq_int_ref = which(seq_int > -1)

  rslt=list(
    seq_ungapped=seq[seq_int > -1],
    seq_int_ungapped=seq_int_ungapped,
    seq_int_ref=seq_int_ref
  )
  return(rslt)
}



bench_seqid=function(seq,seqref){
  return(sum(seq==seqref)/length(seqref))
}


bench_pair=function(seq,seqref,ctref, debug=FALSE){

  npair=sum(ctref$j>0)

  pairs=paste0(seq[ctref$i[ctref$j>0]],seq[ctref$j[ctref$j>0]])

  pairs=toupper(pairs)

  if(debug){
    print(paste(pairs))
  }

  return(sum(pairs %in% RNASSP::energy2)/npair)
}

bench_aln=function(seq,seqref,ctref,debug=FALSE){
  seqid=bench_seqid(seq,seqref)
  pairid=bench_pair(seq,seqref,ctref,debug)
  return(c(
    seqid=seqid,
    pairid=pairid
  ))
}

a2b2seq= function(a2b_1b,seq,mrf_len,type=c("c","s")){

  type=match.arg(type)

  a2b=a2b_1b
  seq_aln = character(mrf_len)
  seq_aln[] = "-"

  seq_aln[a2b[a2b>0]] = seq[a2b > 0]

  if (type=="s"){
    seq_aln = paste0(seq_aln,collapse = "")
  }

  return(seq_aln)

}

#' align sequence to mrf
#'
#' @param seq sequence
#' @param mrf mrf read by read_mrf
#' @param iteration number of iterations
#' @param wt_h weight of field term H
#' @param wt_j weight of coupling term J
#' @param debug verbose
#'
#' @return a2b, mapping index to sequence. Length is the same as mrf
#' @export
#'
align_seq2mrf = function(seq, mrf,iteration=20,wt_h=1.0,wt_j=1.0,debug=TRUE) {

  exp_seq = encode_seq(seq)

  SCO_init = ini_SCO_simple(exp_seq$seq_int_ungapped,
                            mrf_h = mrf$mrf_h)
  SCO_mod = mod_SCO(
    SCO_init,
    iteration = iteration,
    exp_seq$seq_int_ungapped,
    mrf_mat = mrf$mrf_mat,
    mrf_h = mrf$mrf_h,
    wt_h = wt_h,
    wt_j = wt_j,
    DEBUG = debug
  )

  a2b=align(SCO_mod,gap_ext = 0.1,gap_open = -1)
  return(a2b)
}

align_seq2mrf_mtx = function(seq, mrf_mat, mrf_h,iteration=20,wt_h=1.0,wt_j=1.0,debug=TRUE) {

  exp_seq = encode_seq(seq)

  SCO_init = ini_SCO_simple(exp_seq$seq_int_ungapped,
                            mrf_h = mrf$mrf_h)
  SCO_mod = mod_SCO(
    SCO_init,
    iteration = iteration,
    exp_seq$seq_int_ungapped,
    mrf_mat = mrf_mat,
    mrf_h = mrf_h,
    wt_h = wt_h,
    wt_j = wt_j,
    DEBUG = debug
  )

  a2b=align(SCO_mod,gap_ext = 0.1,gap_open = -1)
  return(a2b)
}

#' bencha2b
#'
#' @param a2b_0b a2b 0 based
#' @param seq sequence
#' @param mrf mrf from read_mrf
#' @param seq_ref reference sequence
#' @param ct_ref reference ct with BP
#'
#' @return named vector of benchmark metrics
#' @export
#'
bench_a2b = function(a2b_0b, # a2b output of Rcpp, 0 based
                     seq,
                     mrf,
                     seq_ref = NULL,
                     ct_ref = NULL) {
  exp_seq = encode_seq(seq)
  rslt_aln = score_aln(a2b_0b, exp_seq$seq_int_ungapped, mrf$mrf_mat, mrf$mrf_h, DEBUG = FALSE)


  rslt = c(mrf_total = rslt_aln[1]+rslt_aln[2],
          mrf_single = rslt_aln[1],
           mrf_pair = rslt_aln[2])

  a2b_1b=a2b_0b+1

  seq_final = a2b2seq(a2b_1b, exp_seq$seq_ungapped, mrf$len)

  if (!is.null(seq_ref)) {
    benchaln = bench_aln(seq_final, seqref = seqs[[1]], ctref = ct_ref)
    rslt=c(rslt,benchaln)
  }

  return(rslt)


}

#' sto2ss read SS_cons from stockholm format MSA
#'
#' If multiple lines, paste them together.
#'
#' @param fsto file path to stockholm file with GC SS_cons field
#'
#' @return ss2 Secondary structure in dot bracket format
#' @export
#'
sto2ss=function(fsto){
  ss=system2("cat",args=paste0(fsto,"|grep '^#=GC SS_cons' |grep -v SS_cons2 |awk '{print $3}' "),stdout = TRUE)
  ss2=paste0(ss,collapse = "")
  return(ss2)
}

#' ss2pairs
#'
#' Convert ss to a dataframe with pair info.
#'
#' @param ss Character, Secondary structure in dot bracket format
#'
#' @return dataframe
#'            - id: 1:length(ss)
#'            - pair: paired id, otherwise 0
#'            - ss: character
#' @export
#'
ss2pairs = function(ss) {
  if(length(ss)==1){
    ss=seqinr::s2c(ss)
  }

  pairs = integer(length = length(ss))
  stack1 = integer() # <>
  stack2 = integer() # {}
  stack3 = integer() # ()
  stack4 = integer() # []
  stack5 = integer() # Aa
  stack6 = integer() # Bb
  stack7 = integer() # Cc
  stack8 = integer() # Dd
  for (i in 1:length(ss)) {
    if (ss[i] == "<") {
      stack1 = append(stack1, i)
    }
    else if (ss[i] == ">") {
      pairs[tail(stack1, 1)] = i
      stack1 = head(stack1, -1)
    }
    else if (ss[i] == "A") {
      stack2 = append(stack2, i)
    }
    else if (ss[i] == "a") {
      pairs[tail(stack2, 1)] = i
      stack2 = head(stack2, -1)
    }
    else if (ss[i] == "B") {
      stack3 = append(stack3, i)
    }
    else if (ss[i] == "b") {
      pairs[tail(stack3, 1)] = i
      stack3 = head(stack3, -1)
    }
    else if (ss[i] == "C") {
      stack4 = append(stack4, i)
    }
    else if (ss[i] == "c") {
      pairs[tail(stack4, 1)] = i
      stack4 = head(stack4, -1)
    }
    else if (ss[i] == "D") {
      stack5 = append(stack5, i)
    }
    else if (ss[i] == "d") {
      pairs[tail(stack5, 1)] = i
      stack5 = head(stack5, -1)
    }
    else if (ss[i] == "{") {
      stack6 = append(stack6, i)
    }
    else if (ss[i] == "}") {
      pairs[tail(stack6, 1)] = i
      stack6 = head(stack6, -1)
    }
    else if (ss[i] == "(") {
      stack7 = append(stack7, i)
    }
    else if (ss[i] == ")") {
      pairs[tail(stack7, 1)] = i
      stack7 = head(stack7, -1)
    }
    else if (ss[i] == "[") {
      stack8 = append(stack8, i)
    }
    else if (ss[i] == "]") {
      pairs[tail(stack8, 1)] = i
      stack8 = head(stack8, -1)
    }
  }

  pairidx=which(pairs>0)

  pairs[pairs[pairidx]]=pairidx

  df=data.frame(
    id=1:length(ss),
    pair=pairs,
    ss=ss
  )

  return(df)
}


#' Convert a2m formatted multiple sequence alignment to seqidx, a reference-to-sequence index
#'
#' In an a2m formatted MSA, lowercase residues are insertions, - are gaps, uppercase residues are aligned.
#'
#' @param seq a sequence in a MSA, read from seqinr::read.fasta
#' @param idx_aln integer vector, reference columns. if NULL, it's guessed from the uppercase alignment.
#'
#' @return seqidx, integer vector, pointer from reference columns to sequence.
#' e.g. seqidx[i]=5 means 5th residue is aligned to i.
#' 0 means a gap.
#'
#'
msa_a2m2seqidx=function(seq, idx_aln=NULL){
  idx_seq=which(seq %in% c("A","U","C","G","N","R","M","a","u","c","g","n","r")) # seq -> full alignment

  if (is.null(idx_aln)){
    idx_aln=which(seq %in% c("A","U","C","G","N","R","M","-")) # ref->full alignment
  }

  idx_seqtmp=numeric(length(seq))
  idx_seqtmp[idx_seq]=1:length(idx_seq) # encode seq pointer, full alignment-> seq

  seqidx=idx_seqtmp[idx_aln] # extract aligned seq pointer, ref -> seq
  return(seqidx)
}

msa_a2m2seqidx_all = function(seqs, idx_aln=NULL){
  idx_seqidx=lapply(seqs, function(aseq){
    seq_idx = msa_a2m2seqidx(aseq, idx_aln=idx_aln)
    return(seq_idx)
  })

  seqidx_mat= do.call(rbind,idx_seqidx)
}

get_idx_select = function(dfref) {
  idx_select =
    list(
      col_PK = dfref$id_ref[dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")],
      col_pair = dfref$id_ref[dfref$pair > 0 &
                                (!dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d"))],
      col_loop = dfref$id_ref[dfref$pair == 0],
      col_all = dfref$id_ref
    )
}

read_dfref= function(fsto) {

  ss = RNAmrf::sto2ss(fsto) # seedss
  refs=RNAmrf:::sto2ref(fsto)
  refs.c=seqinr::s2c(refs)
  idx_ref = which(refs.c != ".")

  df = RNAmrf::ss2pairs(ss) # seed ss

  # idx_ref = which(df$ss != ".")
  df$id_ref = 0
  df$id_ref[idx_ref] = 1:length(idx_ref)


  df$id_ref_pair=df$id_ref[match(df$pair,df$id)]
  df$id_ref_pair[is.na(df$id_ref_pair)]=0
  # df$id_mrf=df$id_ref

  dfref = df[df$id_ref > 0, ]
  return(dfref)
}
sto2ref=function(fsto){
  refs=system2("cat",args=paste0(fsto,"|grep '^#=GC RF' |awk '{print $3}' "),stdout = TRUE)
  refs2=paste0(refs,collapse = "")
  return(refs2)
}

msa_a2m2seqidx_mrf=function(seq, dfref){
  seqidx=RNAmrf:::msa_a2m2seqidx(seq)
  seqidx_ref=numeric(length=length(dfref$id_ref))
  seqidx_ref[dfref$id_mrf>0]=seqidx[dfref$id_mrf[(dfref$id_mrf>0)]]
  return(seqidx_ref)
}

a2m2a2b=function(seq, idx_aln=NULL){
  idx_seq=which(seq %in% c("A","U","C","G","N","R","M","a","u","c","g","n","r")) # seq -> full alignment
  if (is.null(idx_aln)){
    idx_aln=which(seq %in% c("A","U","C","G","N","R","M","-")) # ref->full alignment
  }

  idx_seqtmp=numeric(length(seq))
  idx_seqtmp[idx_aln]=1:length(idx_aln)
  a2b=idx_seqtmp[idx_seq]
  return(a2b)
}

bench_pair = function(seqidx.test, seqidx.ref, pair){

  idx_A=which(pair>0)
  idx_B=pair[which(pair>0)]

  true_pair = sum(seqidx.ref[,idx_A]>0)


  predicted_pair=sum(seqidx.test[,idx_A] == seqidx.ref[,idx_A] &
                       seqidx.test[,idx_B] == seqidx.ref[,idx_B] &
                       seqidx.ref[,idx_A] >0)

  sens=predicted_pair/true_pair
  return(sens)
}

bench_by_range = function(seqidx.test, seqidx.ref, range=NULL){

  if(is.null(range)) range=1:ncol(seqidx.ref)

  seqidx.test.r=seqidx.test[,range]
  seqidx.ref.r=seqidx.ref[,range]

  sum(seqidx.ref.r == seqidx.test.r & seqidx.ref.r > 0)/length(which(seqidx.ref.r > 0))
}

bench_dfref= function(seqidx.test, seqidx.ref, dfref){
  idx_select = get_idx_select(dfref)

  rslt=list()

  pair_all=dfref$id_ref_pair
  pair_pk=pair_all
  pair_pk[!dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")]=0

  pair_nonpk=pair_all
  pair_nonpk[dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")]=0


  rslt$col_all=bench_by_range(seqidx.test,seqidx.ref)
  rslt$col_loop=bench_by_range(seqidx.test,seqidx.ref,idx_select$col_loop)
  rslt$pair_all=bench_pair(seqidx.test,seqidx.ref, pair_all)
  rslt$pair_pk=bench_pair(seqidx.test,seqidx.ref, pair_pk)
  rslt$pair_nonpk=bench_pair(seqidx.test,seqidx.ref, pair_nonpk)

  return(unlist(rslt))
}

# get_idx_select = function(dfref) {
#   idx_select =
#     list(
#       col_PK = dfref$id_ref[dfref$ss %in% c("A", "a") & (dfref$id_mrf > 0)],
#       col_pair = dfref$id_ref[dfref$pair > 0 &
#                                 (!dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")) &
#                                 (dfref$id_mrf > 0)],
#       col_loop = dfref$id_ref[dfref$pair == 0 &
#                                 (dfref$id_mrf > 0)],
#       col_all = dfref$id_ref[dfref$id_mrf > 0]
#     )
# }


bench_idx_select = function(seqidx.test,seqidx.ref,idx_select){
  sapply(idx_select, function(a_range) {
    bench_by_range(seqidx.test,seqidx.ref,a_range)
  })
}

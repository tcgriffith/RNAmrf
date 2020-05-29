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

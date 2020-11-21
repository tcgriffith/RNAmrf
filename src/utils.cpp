#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double str_sim(CharacterVector v1,CharacterVector v2){
  double len =v1.length();
  double sim=sum(v1==v2)/len;
  return(sim);
}

// [[Rcpp::export]]
NumericVector calc_wt(CharacterMatrix seqmat,double cutoff=0.8){
  int nseq=seqmat.nrow();
  NumericVector wt(nseq);

  for (int i=0;i < seqmat.nrow();i++){
    CharacterVector vi = seqmat(i,_);
    for(int j=i+1;j<seqmat.nrow();j++){
      CharacterVector vj = seqmat(j,_);
      double sim = str_sim(vi,vj);
      if (sim > cutoff){
        wt[i]++;
        wt[j]++;
      }
    }
  }

  return(wt+1);
}

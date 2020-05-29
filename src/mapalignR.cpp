#include <Rcpp.h>
// #include <math.h>
// #include <algorithm>
#include <iostream>
using namespace Rcpp;

// bool DEBUG=true;

double exp_fast(double x){
  // WARNING fails if |x| > 1024
  //https://codingforspeed.com/using-faster-exponential-approximation/
  x = 1 + x/1024;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}
// [[Rcpp::export]]
double gaussian(double mean, double stdev, double x){return exp_fast(-pow((x - mean),2)/(2*(pow(stdev,2))));}


// [[Rcpp::export]]

NumericVector test(NumericVector x){
	int n=x.size();
	NumericVector out(n);

	out[0]= x[0];
	for (int i = 0; i < n; ++i)
	{
		out[i]=out[i-1]+x[i];
	}

	return out;
}


// [[Rcpp::export]]

double Falign(NumericMatrix sco_mtx){

	int rows=sco_mtx.nrow();
	int cols=sco_mtx.ncol();



    double max_sco = 0;
    NumericMatrix sco(rows+1,cols+1);
    for (int i = 1; i <= rows; i++){
        for (int j = 1; j <= cols; j++){
            double A = sco(i-1,j-1) + sco_mtx(i-1,j-1);
            double D = sco(i-1,j);
            double R = sco(i,j-1);

            if(A >= R){if(A >= D){sco(i,j) = A;}else{sco(i,j) = D;}}
            else{if(R >= D){sco(i,j) = R;}else{sco(i,j) = D;}}

            if(sco(i,j) > max_sco){max_sco = sco(i,j);}
        }
    }
    return(max_sco);
}

int id2_to_id1(int i, int j, int dim_j){
    int id1=(j-1) *dim_j + i;
    return(id1);
}


IntegerVector id1_to_id2(int id1, int dim_i){
    IntegerVector id2(2);

    int mod=id1 % dim_i;

    if (mod == 0){
      id2(0) = dim_i;
      id2(1) = id1 / dim_i;
    } else{
      id2(0) = mod;
      id2(1) = id1 / dim_i +1;
    }

    return(id2);

}



// [[Rcpp::export]]

double retrieve_matj(int i, int a, int j, int b, NumericMatrix mat_j, int len_a){
    return(mat_j( (i)*len_a + a, (j)*len_a + b));
}

double retrieve_mat_h(int i, int a, NumericMatrix mat_h){
  return(mat_h(i,a));
}


// double sepw(double sep){if(sep <= 3){return 0.50;}else if(sep == ){return 0.75;}else{return 1.00;}}
double sepw(double sep){if(sep <= 3){return 0.50;}else{return 1.00;}}

// [[Rcpp::export]]

NumericMatrix ini_SCO(IntegerVector seq,
                      NumericMatrix mrf_mat,
                      NumericMatrix mrf_h,
                      bool DEBUG=false) {


  // SCO = matrix(0, nrow = length(seq), ncol = mrf_len)
  int mrf_len = mrf_h.nrow();
  double sep_x=2;
  double sep_y=1;
  NumericMatrix SCO(seq.length(), mrf_len);

  for (int ai = 0; ai < seq.length(); ++ai){
    int nt_ai=seq(ai);

    for (int bi = 0; bi < mrf_len; ++bi){

      NumericMatrix M(seq.length(), mrf_len);
        for (int aj = 0; aj < seq.length(); aj++){
            for (int bj=0;bj < mrf_len;bj++){

              if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // same direction
                int nt_aj=seq(aj);
                double score_a2b;

                if (bi > bj){
                    score_a2b=retrieve_matj(bj,nt_aj, bi,nt_ai,mrf_mat, 5)+mrf_h(bi,nt_ai)+mrf_h(bj,nt_aj);
                }
                else{
                    score_a2b=retrieve_matj(bi,nt_ai,bj,nt_aj,mrf_mat,5)+mrf_h(bi,nt_ai)+mrf_h(bj,nt_aj); // mrf_mat is symmetric
                }



                int sep_a = abs(ai - aj);
                int sep_b = abs(bi - bj);
                int sep_D = abs(sep_a - sep_b);
                int sep_M = std::min(sep_a, sep_b);
                double sep_std = sep_y * (1 + pow(sep_M - 2,sep_x));

                if (DEBUG){
                  char out[100];
                  std::sprintf(out,"%d %d %f", ai,bi, sep_std);

                  Rcpp::Rcerr << out << std::endl;

                }

                if (sep_D/sep_std <6){
                    M(aj,bj) =  score_a2b *  sepw(sep_M) * gaussian(0,sep_std,sep_D);
                }else{
                    M(aj,bj) = 0;
                }

              }
            }
        }

        SCO(ai,bi)=Falign(M);
    }
  }

  return(SCO);

}

// [[Rcpp::export]]
NumericMatrix ini_SCO_simple(IntegerVector seq, NumericMatrix mrf_h){

  int mrf_len = mrf_h.nrow();
  NumericMatrix SCO(seq.length(), mrf_len);

  // NumericMatrix SCO_h(seq.length(),mrf_len);

  for (int ai = 0; ai < seq.length(); ++ai){
    for (int bi = 0; bi < mrf_len; ++bi){
      int nt_ai=seq(ai);
      SCO(ai,bi) = mrf_h(bi, nt_ai);
    }
  }

  return(SCO);

}


// [[Rcpp::export]]
IntegerVector align(NumericMatrix sco_mtx, double gap_ext,double gap_open){
  // LOCAL_ALIGN
  // Start    0
  // [A]lign  1
  // [D]own   2
  // [R]ight  3

  double max_sco = 0;
  int rows = sco_mtx.nrow();
  int cols = sco_mtx.ncol();

  // vec_int a2b(rows,-1);
  IntegerVector a2b(rows,-1);

  NumericMatrix sco(rows+1,cols+1);
  IntegerMatrix label(rows+1,cols+1);

  int max_i = 0;int max_j = 0;
  for (int i = 1; i <= rows; i++){

    for (int j = 1; j <= cols; j++){
      double A = sco(i-1,j-1) + sco_mtx(i-1,j-1);
      double D = sco(i-1,j);
      double R = sco(i,j-1);

      if(label(i-1,j) == 1){D += gap_open;}else{D += gap_open * gap_ext;}
      if(label(i,j-1) == 1){R += gap_open;}else{R += gap_open * gap_ext;}

      if(A <= 0 and D <= 0 and R <= 0){label(i,j) = 0;sco(i,j) = 0;}
      else{
        if(A >= R){if(A >= D){label(i,j) = 1;sco(i,j) = A;}else{label(i,j) = 2;sco(i,j) = D;}}
        else{if(R >= D){label(i,j) = 3;sco(i,j) = R;}else{label(i,j) = 2;sco(i,j) = D;}}
        if(sco(i,j) > max_sco){max_i = i;max_j = j;max_sco = sco(i,j);}
      }
    }
  }
  int i = max_i;int j = max_j;
  while(1){
    if(label(i,j) == 0){break;}
    else if(label(i,j) == 1){a2b[i-1] = j-1;i--;j--;}
    else if(label(i,j) == 2){i--;}
    else if(label(i,j) == 3){j--;}
  }
  // std::cerr << "aln score: " << max_sco <<std::endl;
  return(a2b);
}


// [[Rcpp::export]]
NumericMatrix align_C_mat(NumericMatrix sco_mtx, double gap_ext,double gap_open){
  // LOCAL_ALIGN
  // Start    0
  // [A]lign  1
  // [D]own   2
  // [R]ight  3

  double max_sco = 0;
  int rows = sco_mtx.nrow();
  int cols = sco_mtx.ncol();

  // vec_int a2b(rows,-1);
  IntegerVector a2b(rows,-1);

  NumericMatrix sco(rows+1,cols+1);
  IntegerMatrix label(rows+1,cols+1);

  int max_i = 0;int max_j = 0;
  for (int i = 1; i <= rows; i++){

    for (int j = 1; j <= cols; j++){
      double A = sco(i-1,j-1) + sco_mtx(i-1,j-1);
      double D = sco(i-1,j);
      double R = sco(i,j-1);

      if(label(i-1,j) == 1){D += gap_open;}else{D += gap_open * gap_ext;}
      if(label(i,j-1) == 1){R += gap_open;}else{R += gap_open * gap_ext;}

      if(A <= 0 and D <= 0 and R <= 0){label(i,j) = 0;sco(i,j) = 0;}
      else{
        if(A >= R){if(A >= D){label(i,j) = 1;sco(i,j) = A;}else{label(i,j) = 2;sco(i,j) = D;}}
        else{if(R >= D){label(i,j) = 3;sco(i,j) = R;}else{label(i,j) = 2;sco(i,j) = D;}}
        if(sco(i,j) > max_sco){max_i = i;max_j = j;max_sco = sco(i,j);}
      }
    }
  }
  int i = max_i;int j = max_j;
  while(1){
    if(label(i,j) == 0){break;}
    else if(label(i,j) == 1){a2b[i-1] = j-1;i--;j--;}
    else if(label(i,j) == 2){i--;}
    else if(label(i,j) == 3){j--;}
  }
  // std::cerr << "aln score: " << max_sco <<std::endl;
  return(sco);
}


// [[Rcpp::export]]
IntegerMatrix align_C_lab(NumericMatrix sco_mtx, double gap_ext,double gap_open){
  // LOCAL_ALIGN
  // Start    0
  // [A]lign  1
  // [D]own   2
  // [R]ight  3

  double max_sco = 0;
  int rows = sco_mtx.nrow();
  int cols = sco_mtx.ncol();

  // vec_int a2b(rows,-1);
  IntegerVector a2b(rows,-1);

  NumericMatrix sco(rows+1,cols+1);
  IntegerMatrix label(rows+1,cols+1);

  int max_i = 0;int max_j = 0;
  for (int i = 1; i <= rows; i++){

    for (int j = 1; j <= cols; j++){
      double A = sco(i-1,j-1) + sco_mtx(i-1,j-1);
      double D = sco(i-1,j);
      double R = sco(i,j-1);

      if(label(i-1,j) == 1){D += gap_open;}else{D += gap_open * gap_ext;}
      if(label(i,j-1) == 1){R += gap_open;}else{R += gap_open * gap_ext;}

      if(A <= 0 and D <= 0 and R <= 0){label(i,j) = 0;sco(i,j) = 0;}
      else{
        if(A >= R){if(A >= D){label(i,j) = 1;sco(i,j) = A;}else{label(i,j) = 2;sco(i,j) = D;}}
        else{if(R >= D){label(i,j) = 3;sco(i,j) = R;}else{label(i,j) = 2;sco(i,j) = D;}}
        if(sco(i,j) > max_sco){max_i = i;max_j = j;max_sco = sco(i,j);}
      }
    }
  }
  int i = max_i;int j = max_j;
  while(1){
    if(label(i,j) == 0){break;}
    else if(label(i,j) == 1){a2b[i-1] = j-1;i--;j--;}
    else if(label(i,j) == 2){i--;}
    else if(label(i,j) == 3){j--;}
  }
  // std::cerr << "aln score: " << max_sco <<std::endl;
  return(label);
}


// [[Rcpp::export]]
IntegerVector align_PSgap(NumericMatrix sco_mtx, double gap_ext,NumericVector gap_open){
  // LOCAL_ALIGN
  // Start    0
  // [A]lign  1
  // [D]own   2
  // [R]ight  3

  double max_sco = 0;
  int rows = sco_mtx.nrow();
  int cols = sco_mtx.ncol();

  // vec_int a2b(rows,-1);
  IntegerVector a2b(rows,-1);

  NumericMatrix sco(rows+1,cols+1);
  IntegerMatrix label(rows+1,cols+1);

  int max_i = 0;int max_j = 0;
  for (int i = 1; i <= rows; i++){

    for (int j = 1; j <= cols; j++){
      double A = sco(i-1,j-1) + sco_mtx(i-1,j-1);
      double D = sco(i-1,j);
      double R = sco(i,j-1);

      if(label(i-1,j) == 1){D += gap_open[i];}else{D += gap_open[i] * gap_ext;}
      if(label(i,j-1) == 1){R += gap_open[i];}else{R += gap_open[i] * gap_ext;}

      if(A <= 0 and D <= 0 and R <= 0){label(i,j) = 0;sco(i,j) = 0;}
      else{
        if(A >= R){if(A >= D){label(i,j) = 1;sco(i,j) = A;}else{label(i,j) = 2;sco(i,j) = D;}}
        else{if(R >= D){label(i,j) = 3;sco(i,j) = R;}else{label(i,j) = 2;sco(i,j) = D;}}
        if(sco(i,j) > max_sco){max_i = i;max_j = j;max_sco = sco(i,j);}
      }
    }
  }
  int i = max_i;int j = max_j;
  while(1){
    if(label(i,j) == 0){break;}
    else if(label(i,j) == 1){a2b[i-1] = j-1;i--;j--;}
    else if(label(i,j) == 2){i--;}
    else if(label(i,j) == 3){j--;}
  }
  // std::cerr << "aln score: " << max_sco <<std::endl;
  return(a2b);
}


// [[Rcpp::export]]
IntegerVector align_PSgap2(NumericMatrix sco_mtx, double gap_ext,NumericVector gap_ins, NumericVector gap_del){
  // LOCAL_ALIGN
  // Start    0
  // [A]lign  1
  // [D]own   2 insert
  // [R]ight  3 deletion

  double max_sco = 0;
  int rows = sco_mtx.nrow();
  int cols = sco_mtx.ncol();

  // vec_int a2b(rows,-1);
  IntegerVector a2b(rows,-1);

  NumericMatrix sco(rows+1,cols+1);
  IntegerMatrix label(rows+1,cols+1);

  int max_i = 0;int max_j = 0;
  for (int i = 1; i <= rows; i++){

    for (int j = 1; j <= cols; j++){
      double A = sco(i-1,j-1) + sco_mtx(i-1,j-1);
      double D = sco(i-1,j);
      double R = sco(i,j-1);

      if(label(i-1,j) == 1){D += gap_ins[i];}else{D += gap_ins[i] * gap_ext;}
      if(label(i,j-1) == 1){R += gap_del[i];}else{R += gap_del[i] * gap_ext;}

      if(A <= 0 and D <= 0 and R <= 0){label(i,j) = 0;sco(i,j) = 0;}
      else{
        if(A >= R){if(A >= D){label(i,j) = 1;sco(i,j) = A;}else{label(i,j) = 2;sco(i,j) = D;}}
        else{if(R >= D){label(i,j) = 3;sco(i,j) = R;}else{label(i,j) = 2;sco(i,j) = D;}}
        if(sco(i,j) > max_sco){max_i = i;max_j = j;max_sco = sco(i,j);}
      }
    }
  }
  int i = max_i;int j = max_j;
  while(1){
    if(label(i,j) == 0){break;}
    else if(label(i,j) == 1){a2b[i-1] = j-1;i--;j--;}
    else if(label(i,j) == 2){i--;}
    else if(label(i,j) == 3){j--;}
  }
  // std::cerr << "aln score: " << max_sco <<std::endl;
  return(a2b);
}




//[[Rcpp::export]]
DoubleVector score_aln(IntegerVector a2b,IntegerVector seq,  NumericMatrix mrf_mat,NumericMatrix mrf_h, bool DEBUG=false){
  DoubleVector score_hj(2);

  // int mrf_len=mrf_h.nrow();

  double score_contact = 0;
  double score_single = 0;

  for(int ai=0; ai <seq.length(); ai++){
    int ai2bi= a2b[ai];
    int nt_ai =seq(ai);
    if (ai2bi >-1){
      score_single += mrf_h(ai2bi, nt_ai);

      for (int aj=0;aj <seq.length();aj++){

        if (aj<ai) continue;

        int aj2bj = a2b[aj];
        if (aj2bj > -1){
          int nt_aj=seq(aj);

          double score_a2b;
          if (ai2bi > aj2bj){
            score_a2b=retrieve_matj(aj2bj,nt_aj, ai2bi,nt_ai,mrf_mat, 5);
          }
          else{
            score_a2b=retrieve_matj(ai2bi,nt_ai,aj2bj,nt_aj,mrf_mat, 5);
          }
          score_contact +=score_a2b;
        }

      }
    }
  }

  score_hj(0) = score_single;
  score_hj(1) = score_contact;

  if (DEBUG){
    char out[100];
    std::sprintf(out,"total: %.2f single: %.2f contact: %.2f", score_single+score_contact, score_single,score_contact);
    Rcpp::Rcerr << out << std::endl;
  }


  return(score_hj);
}



// [[Rcpp::export]]
NumericMatrix mod_SCO(NumericMatrix SCO,
                      int iteration,
                      IntegerVector seq,
                      NumericMatrix mrf_mat,
                      NumericMatrix mrf_h,
                      double wt_h,
                      double wt_j,
                      double gap_o,
                      double gap_e,
                      bool DEBUG=false){

      // iterate
    IntegerVector a2b_tmp;

    NumericMatrix SCO_cln = clone(SCO);

    for(int it=0; it < iteration; it++)
    {
        // align
        a2b_tmp = align(SCO_cln,gap_e,gap_o);

        score_aln(a2b_tmp, seq, mrf_mat, mrf_h,DEBUG);

        // update similarity matrix
        double IT = (double)it + 1;
        double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
        for(int ai=0; ai < SCO_cln.nrow(); ai++){ // go through columns (vec_a) in map_a that has contacts
            // int ai = vec_a[a];
            int nt_ai=seq(ai);
            for(int bi=0; bi < SCO_cln.ncol(); bi++){ // go through columns (vec_b) in map_b that has contacts
                // int bi = vec_b[b];
                double sco_contact = 0;
                double sco_single = 0;
                for(int aj=0; aj < SCO_cln.nrow(); aj++){ // go through contacts in vec_a

                    if (aj == ai) continue;
                    // int aj = vec_a_i[ai,n];
                    int bj = a2b_tmp[aj]; // get mapping
                    int nt_aj=seq(aj);

                    if(bj != -1){ // if mapping exists
                        if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
                            double sep_M = std::min(abs(ai-aj),abs(bi-bj));
                            // sco_contact += mtx_a[ai,aj] * mtx_b[bi,bj] * sepw(sep_M);

                            // int nt_aj=seq(aj);
                            double score_a2b;

                            if (bi > bj){
                                score_a2b=retrieve_matj(bj,nt_aj, bi,nt_ai,mrf_mat, 5);
                            }
                            else{
                                score_a2b=retrieve_matj(bi,nt_ai,bj,nt_aj,mrf_mat, 5);
                            }

                            sco_contact = sco_contact + score_a2b * sepw(sep_M);


                            // sco_contact = score_a2b *  sepw(sep_M);
                            // sco_contact = 0;
                        }
                    }
                }
                if (ai==2 ) {
                  // std::cerr << bi << " " << sco_contact << std::endl;
                }

                // double wt_single = 0.0;
                // double wt_contact = 1.0;
                sco_single = mrf_h(bi, nt_ai);
                // SCO_cln(ai,bi) = s1 *SCO_cln(ai,bi) + s2 *sco_contact;
                SCO_cln(ai,bi) = s1*SCO_cln(ai,bi) + s2*(sco_contact*wt_j+sco_single*wt_h);

            }
        }
    }

  return(SCO_cln);

}

// [[Rcpp::export]]
NumericMatrix mod_SCO_PSgap(NumericMatrix SCO,
                      int iteration,
                      IntegerVector seq,
                      NumericMatrix mrf_mat,
                      NumericMatrix mrf_h,
                      double wt_h,
                      double wt_j,
                      NumericVector gap_o, // position specific gaps
                      double gap_e,
                      bool DEBUG=false){

  // iterate
  IntegerVector a2b_tmp;

  NumericMatrix SCO_cln = clone(SCO);

  for(int it=0; it < iteration; it++)
  {
    // align
    a2b_tmp = align_PSgap(SCO_cln,gap_e,gap_o);

    score_aln(a2b_tmp, seq, mrf_mat, mrf_h,DEBUG);

    // update similarity matrix
    double IT = (double)it + 1;
    double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
    for(int ai=0; ai < SCO_cln.nrow(); ai++){ // go through columns (vec_a) in map_a that has contacts
      // int ai = vec_a[a];
      int nt_ai=seq(ai);
      for(int bi=0; bi < SCO_cln.ncol(); bi++){ // go through columns (vec_b) in map_b that has contacts
        // int bi = vec_b[b];
        double sco_contact = 0;
        double sco_single = 0;
        for(int aj=0; aj < SCO_cln.nrow(); aj++){ // go through contacts in vec_a

          if (aj == ai) continue;
          // int aj = vec_a_i[ai,n];
          int bj = a2b_tmp[aj]; // get mapping
          int nt_aj=seq(aj);

          if(bj != -1){ // if mapping exists
            if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
              double sep_M = std::min(abs(ai-aj),abs(bi-bj));
              // sco_contact += mtx_a[ai,aj] * mtx_b[bi,bj] * sepw(sep_M);

              // int nt_aj=seq(aj);
              double score_a2b;

              if (bi > bj){
                score_a2b=retrieve_matj(bj,nt_aj, bi,nt_ai,mrf_mat, 5);
              }
              else{
                score_a2b=retrieve_matj(bi,nt_ai,bj,nt_aj,mrf_mat, 5);
              }

              sco_contact = sco_contact + score_a2b *sepw(sep_M);


              // sco_contact = score_a2b *  sepw(sep_M);
              // sco_contact = 0;
            }
          }
        }
        if (ai==1) {
          // std::cerr << bi << " " << sco_contact << std::endl;
        }

        // double wt_single = 0.0;
        // double wt_contact = 1.0;
        sco_single = mrf_h(bi, nt_ai);
        // SCO_cln(ai,bi) = s1 *SCO_cln(ai,bi) + s2 *sco_contact;
        SCO_cln(ai,bi) = s1*SCO_cln(ai,bi) + s2*(sco_contact*wt_j+sco_single*wt_h);
      }
    }
  }

  return(SCO_cln);

}


// [[Rcpp::export]]
NumericMatrix mod_SCO_PSgap2(NumericMatrix SCO,
                      int iteration,
                      IntegerVector seq,
                      NumericMatrix mrf_mat,
                      NumericMatrix mrf_h,
                      double wt_h,
                      double wt_j,
                      NumericVector gap_ins, // position specific gaps, insertion
                      NumericVector gap_del, // position specific gaps, deletion
                      double gap_e,
                      bool DEBUG=false){

  // iterate
  IntegerVector a2b_tmp;

  NumericMatrix SCO_cln = clone(SCO);

  for(int it=0; it < iteration; it++)
  {
    // align
    a2b_tmp = align_PSgap2(SCO_cln,gap_e,gap_ins, gap_del);

    score_aln(a2b_tmp, seq, mrf_mat, mrf_h,DEBUG);

    // update similarity matrix
    double IT = (double)it + 1;
    double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
    for(int ai=0; ai < SCO_cln.nrow(); ai++){ // go through columns (vec_a) in map_a that has contacts
      // int ai = vec_a[a];
      int nt_ai=seq(ai);
      for(int bi=0; bi < SCO_cln.ncol(); bi++){ // go through columns (vec_b) in map_b that has contacts
        // int bi = vec_b[b];
        double sco_contact = 0;
        double sco_single = 0;
        for(int aj=0; aj < SCO_cln.nrow(); aj++){ // go through contacts in vec_a

          if (aj == ai) continue;
          // int aj = vec_a_i[ai,n];
          int bj = a2b_tmp[aj]; // get mapping
          int nt_aj=seq(aj);

          if(bj != -1){ // if mapping exists
            if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
              double sep_M = std::min(abs(ai-aj),abs(bi-bj));
              // sco_contact += mtx_a[ai,aj] * mtx_b[bi,bj] * sepw(sep_M);

              // int nt_aj=seq(aj);
              double score_a2b;

              if (bi > bj){
                score_a2b=retrieve_matj(bj,nt_aj, bi,nt_ai,mrf_mat, 5);
              }
              else{
                score_a2b=retrieve_matj(bi,nt_ai,bj,nt_aj,mrf_mat, 5);
              }

              sco_contact = sco_contact + score_a2b *sepw(sep_M);


              // sco_contact = score_a2b *  sepw(sep_M);
              // sco_contact = 0;
            }
          }
        }
        if (ai==1) {
          // std::cerr << bi << " " << sco_contact << std::endl;
        }

        // double wt_single = 0.0;
        // double wt_contact = 1.0;
        sco_single = mrf_h(bi, nt_ai);
        // SCO_cln(ai,bi) = s1 *SCO_cln(ai,bi) + s2 *sco_contact;
        SCO_cln(ai,bi) = s1*SCO_cln(ai,bi) + s2*(sco_contact*wt_j+sco_single*wt_h);
      }
    }
  }

  return(SCO_cln);

}


String rna_i2c_upper(int NT){
  if(NT==0) return("A");
  else if(NT==1) return("U");
  else if(NT==2) return("C");
  else if(NT==3) return("G");
  else return("X");
}

String rna_i2c_lower(int NT){
  if(NT==0) return("a");
  else if(NT==1) return("u");
  else if(NT==2) return("c");
  else if(NT==3) return("g");
  else return("x");
}

// [[Rcpp::export]]
CharacterVector a2b2a2m(IntegerVector a2b, IntegerVector seq, int mrflen){
  int last_idx=-1;

  IntegerVector a2b0b = clone(a2b);
  int inslen=0;

  for (int i=0; i < a2b0b.length();i++){
    if (a2b0b[i]==-1){
       for(int j=i+1;j<a2b0b.length();j++){
         if (a2b0b[j]>last_idx){
           a2b0b[j]++;
         }
       }
       a2b0b[i]=last_idx +1;
      inslen++;
    }
    last_idx = a2b0b[i];
  }

  CharacterVector seqa2m(mrflen+inslen,"-");

  for(int i=0;i<a2b0b.length();i++){
    if (a2b[i]==-1){
      seqa2m[a2b0b[i]]=rna_i2c_lower(seq[i]);
    }
    else{
      seqa2m[a2b0b[i]]=rna_i2c_upper(seq[i]);
    }
  }
  return(seqa2m);
}




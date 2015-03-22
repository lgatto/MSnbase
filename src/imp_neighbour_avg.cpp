#include <Rcpp.h>
using namespace Rcpp;

bool allZero(NumericVector x, double zero){
  for (int i=0; i < x.size(); i++) {
    if (x[i] != zero) return false;
  }
  return true;
}

NumericVector subset(NumericVector x, int i) {
  // seq_len is 1-based
  NumericVector y = x[seq_len(i)-1];
  return y;
}



// [[Rcpp::export]]
NumericMatrix imp_neighbour_avg(NumericMatrix x, double k) {
  // input matrix is expected to have >= 3 columns
  NumericMatrix ans = clone(x);
  int nr = ans.nrow(), nc = ans.ncol();

  for(int i = 0; i < nr; i++) {
    // first and last values are set to 0 if NA
    if (R_IsNA(ans(i, 0))) ans(i, 0) = k;
    if (R_IsNA(ans(i, nc-1))) ans(i, nc-1) = k;
    
    for(int j = 1; j < (nc-1); j++) {
      if (R_IsNA(ans(i,j))) {
	// if the next value is NA and all previous values are 0
	// then we set to 0	
	if (R_IsNA(ans(i,j+1))) {
	  NumericVector v = subset(ans.row(i), j);
	  if (allZero(v, k)) ans(i,j) = k;
	} else { // next is not NA, set to mean of neighbours
	  ans(i,j) = (ans(i,j-1) + ans(i,j+1))/2;
	}
      }
    }
  }    
  return(ans);
}


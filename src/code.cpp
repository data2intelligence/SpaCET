#include <Rcpp.h>
using namespace Rcpp;

//' Multiply a number by two
//' 
//' @param x A single integer.
//' @export
// [[Rcpp::export]]

int cl(int x) {
  int y=0;
  for (int i = 0; i <= x; i++) {
    y=y+i;
  };
  return y;
}
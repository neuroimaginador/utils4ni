#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
IntegerMatrix confusion_matrix(IntegerVector label1, IntegerVector label2) {

  int mx = max(label1);
  int my = max(label2);

  IntegerMatrix result(mx + 1, my + 1);

  for (int idx = 0; idx < label1.size(); idx++) {

    result(label1[idx], label2[idx])++;

  }

  return result;

}

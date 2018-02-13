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
NumericVector mask_values(NumericVector input, double low, double high, double in_val, double out_val) {

  NumericVector res(input.size());
  res.attr("dim") = input.attr("dim");

  for (int i = 0; i < res.size(); i++) {

    if ((input[i] >= low) & ((input[i] <= high))) {

      res[i] = in_val;

    } else {

      res[i] = out_val;

    }

  }

  return res;

}

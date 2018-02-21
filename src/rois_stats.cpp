#include <Rcpp.h>
#include <memory.h>
#include <math.h>
// #include "pointers.h"
// #include "pointers_ops.h"

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector sum_4d(NumericVector values) {

  IntegerVector dims = values.attr("dim");
  IntegerVector new_dims(3);
  for (int d = 0; d < 3; d++) {

    new_dims[d] = dims[d];

  }
  int n_4d = dims[3];

  NumericVector suma(values.size() / n_4d);

  for (int voxel = 0; voxel < suma.size(); voxel++) {

    for (int d = 0; d < n_4d; d++) {

      suma[voxel] += values[d * suma.size() + voxel];

    }

  }

  suma.attr("dim") = new_dims;

  return suma;

}

// [[Rcpp::export]]
NumericVector sum_by_ROI(IntegerVector labelled, NumericVector values) {

  int N = max(labelled);
  NumericVector out(N, 0.0);

  for (int idx = 0; idx < labelled.size(); idx++) {

    if (labelled[idx] > 0) {

      out[labelled[idx] - 1] += values[idx];

    }

  }

  return out;

}

// [[Rcpp::export]]
NumericVector max_by_ROI(IntegerVector labelled, NumericVector values) {

  int N = max(labelled);

  NumericVector out(N, -INFINITY);

  for (int idx = 0; idx < labelled.size(); idx++) {

    if (labelled[idx] > 0) {

      if (values[idx] > out[labelled[idx] - 1]) {

        out[labelled[idx] - 1] = values[idx];

      }

    }

  }

  return out;
}

// [[Rcpp::export]]
NumericVector min_by_ROI(IntegerVector labelled, NumericVector values) {

  return -max_by_ROI(labelled, -values);

}

// [[Rcpp::export]]
IntegerVector count_by_ROI(IntegerVector labelled) {

  int N = max(labelled);

  IntegerVector out(N, 0.0);

  for (int idx = 0; idx < labelled.size(); idx++) {

    if (labelled[idx] > 0) {

      out[labelled[idx] - 1]++;

    }

  }

  return out;

}

// [[Rcpp::export]]
NumericVector mean_by_ROI(IntegerVector labelled, NumericVector values) {

  int N = max(labelled);

  NumericVector out(N, 0.0);

  NumericVector sum = sum_by_ROI(labelled, values);
  IntegerVector count = count_by_ROI(labelled);

  out = sum / (NumericVector)count;


  return out;

}


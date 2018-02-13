#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
// #include "defuzzify.h"
// #include "regularize.h"

// [[Rcpp::plugins("cpp11")]]

#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(p) { if(p) { delete[] (p); (p) = NULL; } }
#endif

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/


void pve_segment(double* data, int nvoxels, double* otsu_estimates, int nclasses, double* segmentation) {

  for (int i = 0; i < nvoxels; i++) {

    if (data[i] <= otsu_estimates[0] ) {

      segmentation[i] = 1;
      continue;

    }

    if (data[i] >= otsu_estimates[nclasses - 1] ) {

      segmentation[(nclasses - 1) * nvoxels + i] = 1;
      continue;

    }

    if (nclasses > 1) {

      // double min_distance = 1.e10;

      for (int k = 0; k < nclasses; k++ ) {
        double m = otsu_estimates[k - 1];
        double M = otsu_estimates[k];

        if (data[i] >= m && data[i] <= M) {

          segmentation[k * nvoxels + i] = (data[i] - m) / (M - m);
          segmentation[(k - 1) * nvoxels + i] = 1 - segmentation[k * nvoxels + i];

        }

        // segmentation[k * nvoxels + i] = fabs(data[i] - otsu_estimates[k]);
        //
        // if (segmentation[k * nvoxels + i] < min_distance) {
        //
        //   min_distance = segmentation[k * nvoxels + i];
        //
        // }

      }

      // for (int k = 0; k < nclasses; k++) {
      //
      //   segmentation[k * nvoxels + i] = exp(- segmentation[k * nvoxels + i] / (0.15 * min_distance + 1.e-15));
      //
      // }

    }


  }

}

void segment(double* data, int ndims, int* dims, int nclasses, double* otsu_estimates, double* segmentation) {

  int nvoxels = 1;
  for (int i = 0; i < ndims; i++) {

    nvoxels *= dims[i];

  }

  pve_segment(data, nvoxels, otsu_estimates, nclasses, segmentation);

}

//[[Rcpp::export]]
NumericVector segmentation(NumericVector image, NumericVector otsu_estimates) {

  IntegerVector dims = image.attr("dim");
  IntegerVector segmentation_dims(dims.size() + 1);

  // int* segmentation_dims = (int*) malloc((dims.size() + 1) * sizeof(int));
  for (int i = 0; i < dims.size(); i++) {

    (segmentation_dims.begin())[i] = dims[i];

  }

  segmentation_dims[dims.size()] = otsu_estimates.size();
  NumericVector segmentation(image.size() * otsu_estimates.size());

  segment(image.begin(), dims.size(), dims.begin(), otsu_estimates.size(), otsu_estimates.begin(), segmentation.begin());
  segmentation.attr("dim") = segmentation_dims;

  return segmentation;

}


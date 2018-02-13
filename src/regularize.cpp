#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

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
//

void regularize(double* input, int ndims, int* dims, double* kernel, int kernel_ndims, int kernel_width, double* output) {

  // We should have (ndims == kernel_ndims) or (ndims == kernel_ndims + 1)
  //

  if (ndims == 3) {

    if (kernel_width % 2 == 1) kernel_width--;
    int search_limit = kernel_width / 2;
    int limit = search_limit + 1;

    for (int x = limit; x < dims[0] - limit; x++) {

      for (int y = limit; y < dims[1] - limit; y++) {

        for (int z = limit; z < dims[2] - limit; z++) {

          int voxel = x + dims[0] * y + dims[0] * dims[1] * z;

          double cumul = 0.0;

          for (int xdisp = -search_limit; xdisp < search_limit; xdisp++) {

            for (int ydisp = -search_limit; ydisp < search_limit; ydisp++) {

              for (int zdisp = -search_limit; zdisp < search_limit; zdisp++) {

                int loc_kernel = (xdisp + search_limit) + kernel_width * (ydisp + search_limit) + kernel_width * kernel_width * (zdisp + search_limit);
                int loc_image = (x + xdisp) + dims[0] * (y + ydisp) + dims[0] * dims[1] * (z + zdisp);
                cumul += input[loc_image] * kernel[loc_kernel];

              }

            }

          }

          output[voxel] = cumul;

        }

      }

    }

  } else {

    if (kernel_width % 2 == 1) kernel_width--;
    int search_limit = kernel_width / 2;
    int limit = search_limit + 1;

    for (int x = limit; x < dims[0] - limit; x++) {

      for (int y = limit; y < dims[1] - limit; y++) {

        for (int z = limit; z < dims[2] - limit; z++) {

          for (int k = 0; k < dims[3]; k++) {

            int voxel = x + dims[0] * y + dims[0] * dims[1] * z +  dims[0] * dims[1] * dims[2] * k;

            double cumul = 0.0;

            for (int xdisp = -search_limit; xdisp < search_limit; xdisp++) {

              for (int ydisp = -search_limit; ydisp < search_limit; ydisp++) {

                for (int zdisp = -search_limit; zdisp < search_limit; zdisp++) {

                  int loc_kernel = (xdisp + search_limit) + kernel_width * (ydisp + search_limit) + kernel_width * kernel_width * (zdisp + search_limit);
                  int loc_image = (x + xdisp) + dims[0] * (y + ydisp) + dims[0] * dims[1] * (z + zdisp) +  dims[0] * dims[1] * dims[2] * k;
                  cumul += input[loc_image] * kernel[loc_kernel];

                }

              }

            }

            output[voxel] = cumul;

          }

        }

      }

    }

  }

}


//[[Rcpp::export]]
NumericVector regularize(NumericVector image, NumericVector kernel) {

  IntegerVector dims = image.attr("dim");
  int ndims = dims.size();

  IntegerVector kernel_dims = kernel.attr("dim");
  int kernel_ndims = kernel_dims.size();
  int kernel_width = kernel_dims[0];

  NumericVector segmentation(image.size());

  regularize(image.begin(), ndims, dims.begin(), kernel.begin(), kernel_ndims, kernel_width, segmentation.begin());

  segmentation.attr("dim") = dims;

  return segmentation;

}

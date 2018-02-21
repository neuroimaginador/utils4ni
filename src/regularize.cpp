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

void regularize(double* input,
                int ndims,
                int* dims,
                double* kernel,
                int kernel_ndims,
                int kernel_width,
                double* output,
                int ncores = 1) {

  // We should have (ndims == kernel_ndims) or (ndims == kernel_ndims + 1)
  //

  if (ndims == 3) {

    if (kernel_width % 2 == 1) kernel_width--;
    int search_limit = kernel_width / 2;
    int limit = search_limit + 1;

    // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores)
    for (int x = 0; x < dims[0]; x++) {

      for (int y = 0; y < dims[1]; y++) {

        for (int z = 0; z < dims[2]; z++) {

          int voxel = x + dims[0] * y + dims[0] * dims[1] * z;

          double cumul = 0.0;
          double cumul_kernel = 0.0;

          for (int xdisp = -search_limit; xdisp < search_limit; xdisp++) {

            for (int ydisp = -search_limit; ydisp < search_limit; ydisp++) {

              for (int zdisp = -search_limit; zdisp < search_limit; zdisp++) {

                int loc_kernel = (xdisp + search_limit) + kernel_width * (ydisp + search_limit) + kernel_width * kernel_width * (zdisp + search_limit);
                int loc_image = (x + xdisp) + dims[0] * (y + ydisp) + dims[0] * dims[1] * (z + zdisp);

                float value = 0;
                if ((x + xdisp >= 0) && (x + xdisp < dims[0]) &&
                    (y + ydisp >= 0) && (y + ydisp < dims[1]) &&
                    (z + zdisp >= 0) && (z + zdisp < dims[2])) {

                  value = input[loc_image] * kernel[loc_kernel];
                  cumul_kernel += kernel[loc_kernel];
                }

                cumul += value;

              }

            }

          }

          output[voxel] = cumul;
          if (cumul_kernel > 0)
            output[voxel] /= cumul_kernel;

        }

      }

    }

  } else {

    if (kernel_width % 2 == 1) kernel_width--;
    int search_limit = kernel_width / 2;
    int limit = search_limit + 1;

    // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores)
    for (int x = 0; x < dims[0]; x++) {

      for (int y = 0; y < dims[1]; y++) {

        for (int z = 0; z < dims[2]; z++) {

          for (int k = 0; k < dims[3]; k++) {

            int voxel = x + dims[0] * y + dims[0] * dims[1] * z +  dims[0] * dims[1] * dims[2] * k;

            double cumul = 0.0;
            double cumul_kernel = 0.0;

            for (int xdisp = -search_limit; xdisp < search_limit; xdisp++) {

              for (int ydisp = -search_limit; ydisp < search_limit; ydisp++) {

                for (int zdisp = -search_limit; zdisp < search_limit; zdisp++) {

                  int loc_kernel = (xdisp + search_limit) + kernel_width * (ydisp + search_limit) + kernel_width * kernel_width * (zdisp + search_limit);
                  int loc_image = (x + xdisp) + dims[0] * (y + ydisp) + dims[0] * dims[1] * (z + zdisp) +  dims[0] * dims[1] * dims[2] * k;

                  float value = 0;
                  if ((x + xdisp >= 0) && (x + xdisp < dims[0]) &&
                      (y + ydisp >= 0) && (y + ydisp < dims[1]) &&
                      (z + zdisp >= 0) && (z + zdisp < dims[2])) {

                    value = input[loc_image] * kernel[loc_kernel];
                    cumul_kernel += kernel[loc_kernel];
                  }

                  cumul += value;

                }

              }

            }

            output[voxel] = cumul;
            if (cumul_kernel > 0)
              output[voxel] /= cumul_kernel;

          }

        }

      }

    }

  }

}


//[[Rcpp::export]]
NumericVector regularize(NumericVector image, NumericVector kernel, int ncores = 1) {

  IntegerVector dims = image.attr("dim");
  int ndims = dims.size();

  IntegerVector kernel_dims = kernel.attr("dim");
  int kernel_ndims = kernel_dims.size();
  int kernel_width = kernel_dims[0];

  NumericVector segmentation(image.size());

  regularize(image.begin(),
             ndims,
             dims.begin(),
             kernel.begin(),
             kernel_ndims,
             kernel_width,
             segmentation.begin(),
             ncores);

  segmentation.attr("dim") = dims;

  return segmentation;

}

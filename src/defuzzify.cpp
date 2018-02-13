#include <Rcpp.h>
// #include "defuzzify.h"
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


void defuzzify(double* fuzzy_voting, int n_voxels_full, int n_labels, int* label_ids, int* segmentation, int ncores = 1) {

  // Compute final segmentation
  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores)
  for (int voxel = 0; voxel < n_voxels_full; voxel++) {

    double max_voting = 0;
    int mylabel = 0;

    for (int lb = 0; lb < n_labels; lb++) {

      if (fuzzy_voting[lb * n_voxels_full + voxel] > max_voting) {

        max_voting = fuzzy_voting[lb * n_voxels_full + voxel];
        mylabel = lb;

      }

    }

    // if (max_voting < 0.5) mylabel = 0;

    segmentation[voxel] = label_ids[mylabel];

  }

}

void defuzzify(double* fuzzy_voting, int n_voxels_full, int n_labels, int* segmentation, int ncores = 1) {

  // Compute final segmentation
  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores)
  for (int voxel = 0; voxel < n_voxels_full; voxel++) {

    double max_voting = 0;
    int mylabel = 0;

    for (int lb = 0; lb < n_labels; lb++) {

      if (fuzzy_voting[lb * n_voxels_full + voxel] > max_voting) {

        max_voting = fuzzy_voting[lb * n_voxels_full + voxel];
        mylabel = lb;

      }

    }

    segmentation[voxel] = mylabel + 1;

  }

}

//[[Rcpp::export]]
IntegerVector defuzzify(NumericVector image, int ncores = 1) {

  IntegerVector dims = image.attr("dim");
  IntegerVector segmentation_dims(dims.size() - 1);

  int nvoxels = 1;

  for (int i = 0; i < dims.size() - 1; i++) {

    segmentation_dims[i] = dims[i];
    nvoxels *= dims[i];

  }

  IntegerVector segmentation(nvoxels);
  int nclasses = (dims.begin())[dims.size() - 1];

  defuzzify(image.begin(), nvoxels, nclasses, segmentation.begin());

  segmentation.attr("dim") = segmentation_dims;

  return segmentation;

}



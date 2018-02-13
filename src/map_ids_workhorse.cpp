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
IntegerVector map_ids_workhorse(IntegerVector x, IntegerVector source, IntegerVector target) {

  int n_voxels = x.length();
  IntegerVector res(n_voxels);

  int n_values = source.length();

  for (int voxel = 0; voxel < n_voxels; voxel++) {

    res[voxel] = 0;

    if (x[voxel] > 0) {

      for (int id = 0; id < n_values; id++) {

        if (x[voxel] == source[id]) {

          res[voxel] = target[id];
          break;

        }

      }

    }

  }

  IntegerVector dims = x.attr("dim");
  res.attr("dim") = dims;

  return(res);

}

// [[Rcpp::export]]
IntegerVector map_extra_classes(IntegerVector x, IntegerVector source, int remaining) {

  int n_voxels = x.length();
  IntegerVector res(n_voxels);

  int n_values = source.length();

  for (int voxel = 0; voxel < n_voxels; voxel++) {

    if (x[voxel] > 0) {

      res[voxel] = remaining;

      for (int id = 0; id < n_values; id++) {

        if (x[voxel] == source[id]) {

          res[voxel] = x[voxel];
          break;

        }

      }

    }

  }

  IntegerVector dims = x.attr("dim");
  res.attr("dim") = dims;

  return(res);

}

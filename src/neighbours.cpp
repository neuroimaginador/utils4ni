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


void get_neighbours_relative_coordinates(IntegerVector array,
                                         int width,
                                         int* neighbours_offsets,
                                         int* locs) {

  IntegerVector dims = array.attr("dim");
  int ndims = dims.size();

  int* extremes = (int*) malloc(ndims * sizeof(int));

  if (width % 2 == 0) width++;

  for (int i = 0; i < ndims; i++) {

    extremes[i] = (width - 1) / 2;

  }

  int neighbourhood_size = 1;
  int* steps = (int*) malloc((ndims + 1) * sizeof(int));
  steps[0] = 1;

  for (int i = 0; i < ndims; i++) {

    neighbourhood_size *= width;
    steps[i + 1] = steps[i] * dims[i];

  }

  // int* locs = (int*) malloc(neighbourhood_size * array->ndims * sizeof(int));

  for (int j = 0; j < neighbourhood_size; j++) {

    if (j == 0) {

      for (int i = 0; i < ndims; i++)
        locs[j + i * neighbourhood_size] = -extremes[i];

    } else {

      locs[j + 0 * neighbourhood_size] = locs[(j - 1) + 0 * neighbourhood_size] + 1;

      for (int i = 0; i < ndims; i++) {

        if (locs[j + i * neighbourhood_size] > extremes[i]) {

          locs[j + i * neighbourhood_size] = -extremes[i];
          locs[j + (i + 1) * neighbourhood_size] = locs[(j - 1) + (i + 1) * neighbourhood_size] + 1;

        } else if (i < (ndims - 1))
          locs[j + (i + 1) * neighbourhood_size] = locs[(j - 1) + (i + 1) * neighbourhood_size];

      }

    }

    neighbours_offsets[j] = 0;

    for (int i = 0; i < ndims; i++)
      neighbours_offsets[j] += locs[j + i * neighbourhood_size] * steps[i];

  }


  free(steps);
  free(extremes);

}

void get_neighbours_indices(int* dims,
                            int ndims,
                            int width,
                            int* neighbours_offsets) {

  int* extremes = (int*) malloc(ndims * sizeof(int));

  if (width % 2 == 0) width++;

  for (int i = 0; i < ndims; i++) {

    extremes[i] = (width - 1) / 2;

  }

  int neighbourhood_size = 1;
  int* steps = (int*) malloc((ndims + 1) * sizeof(int));
  steps[0] = 1;

  for (int i = 0; i < ndims; i++) {

    neighbourhood_size *= width;
    steps[i + 1] = steps[i] * dims[i];

  }

  int* locs = (int*) malloc(neighbourhood_size * ndims * sizeof(int));

  for (int j = 0; j < neighbourhood_size; j++) {

    if (j == 0) {

      for (int i = 0; i < ndims; i++)
        locs[j + i * neighbourhood_size] = -extremes[i];

    } else {

      locs[j + 0 * neighbourhood_size] = locs[(j - 1) + 0 * neighbourhood_size] + 1;

      for (int i = 0; i < ndims; i++) {

        if (locs[j + i * neighbourhood_size] > extremes[i]) {

          locs[j + i * neighbourhood_size] = -extremes[i];
          locs[j + (i + 1) * neighbourhood_size] = locs[(j - 1) + (i + 1) * neighbourhood_size] + 1;

        } else if (i < (ndims - 1))
          locs[j + (i + 1) * neighbourhood_size] = locs[(j - 1) + (i + 1) * neighbourhood_size];

      }

    }

    neighbours_offsets[j] = 0;

    for (int i = 0; i < ndims; i++)
      neighbours_offsets[j] += locs[j + i * neighbourhood_size] * steps[i];

  }

  free(locs);
  free(steps);
  free(extremes);

}


// [[Rcpp::export]]
IntegerVector get_neighbours(NumericVector array, int width) {


  // Rprintf("B\n");
  IntegerVector dims = array.attr("dim");
  int ndims = dims.length();
  int neighbourhood_size = pow(width, ndims);

  IntegerVector output(neighbourhood_size);

  // Rprintf("C\n");
  get_neighbours_indices(dims.begin(), ndims, width, output.begin());

  return output;

}

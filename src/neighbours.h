#include <Rcpp.h>
#include <Rdefines.h>
using namespace Rcpp;

void get_neighbours_relative_coordinates(IntegerVector array,
                                         int width,
                                         int* neighbours_offsets,
                                         int* locs);

void get_neighbours_indices(int* dims,
                            int ndims,
                            int width,
                            int* neighbours_offsets);

IntegerVector get_neighbours(NumericVector array, int width);

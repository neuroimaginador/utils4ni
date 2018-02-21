#include <Rcpp.h>
#include <Rinternals.h>
using namespace Rcpp;

void get_neighbours(int idx, int* neighbours, int* neighs, int n_neighbours);

template <typename T>
void get_values(const T* image, int* idx, T* values, int n);

void get_image_value(const double* image, int* idx, double* values, int n);

void get_label_value(const int* image, int* idx, int* values, int n);

void get_neighbours_relative_coordinates(IntegerVector array,
                                         int width,
                                         int* neighbours_offsets,
                                         int* locs);

void get_neighbours_indices(int* dims,
                            int ndims,
                            int width,
                            int* neighbours_offsets);

IntegerVector get_neighbours(NumericVector array, int width);

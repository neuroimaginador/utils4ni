#include <Rcpp.h>
#include "neighbours.h"
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


int find(int i, int* nodes) {

  while (nodes[i] != i) {

    i = nodes[i];

  }

  return nodes[i];

}


// [[Rcpp::export]]
IntegerVector connected_components(IntegerVector image) {

  IntegerVector labels(image.size(), 0);
  IntegerVector dims = image.attr("dim");
  labels.attr("dim") = dims;
  int ndims = dims.size();

  // First pass
  int nextlabel = 1;
  int* nodes = (int*) malloc((image.size()) / 2 * sizeof(int));
  // memset(labels.begin()->data, 0, labels->n * sizeof(int));

  // Rprintf("B\n");
  int neighbourhood_size = pow(3, ndims);
  int* neighbours_offsets = (int*) malloc(neighbourhood_size * sizeof(int));
  int* neighbours_locs = (int*) malloc(neighbourhood_size * ndims *  sizeof(int));
  int* neighs = (int*) malloc(neighbourhood_size * sizeof(int));
  int* label_neighs = (int*) malloc(neighbourhood_size * sizeof(int));
  int* expanded_idx = (int*) malloc(ndims * sizeof(int));

  // Rprintf("C\n");
  // get_neighbours_indices(image, 3, neighbours_offsets);
  get_neighbours_relative_coordinates(image,
                                      3,
                                      neighbours_offsets,
                                      neighbours_locs);

  // Rprintf("D\n");

  for (int i = 0; i < image.size(); i++) {

    int mylabel = image[i];

    if (mylabel > 0) {

      // double min_label = 1.e10;
      std::vector<int> partial_labels;

      int x_orig = i % dims[0];
      int y_orig = ((i - x_orig) / dims[0]) % dims[1];
      int z_orig = ((i - x_orig) / dims[0] - y_orig) / dims[1];

      // Get actual neighbours
      expanded_idx[0] = x_orig;
      expanded_idx[1] = y_orig;
      expanded_idx[2] = z_orig;
      // expand_index(image, i, expanded_idx);

      for (int j = 0; j < neighbourhood_size; j++) {

        // Check if the jth neighbour is inside

        bool is_inside = true;
        for (int k = 0; k < ndims; k++) {

          int partial_idx = neighbours_locs[j + k * neighbourhood_size] + expanded_idx[k];

          if (partial_idx < 0 || partial_idx >= dims[k]) {

            is_inside = false;
            break;

          }

        }

        if (is_inside &&
            neighbours_offsets[j] != 0 &&
            (image[neighbours_offsets[j] + i] == mylabel) &&
            (labels[neighbours_offsets[j] + i] > 0)) {

          partial_labels.push_back( labels[neighbours_offsets[j] + i]);

        }

      }

      if (partial_labels.size() == 0) {

        nodes[nextlabel] = nextlabel;
        labels[i] = nextlabel;

      } else {

        int ML = 1.e8;

        for (int k = 0; k < (int)partial_labels.size(); k++) {

          if (ML > partial_labels[k]) {

            ML = partial_labels[k];

          }

        }

        labels[i] = (int) ML;
        nodes[nextlabel] = ML;

        while (partial_labels.size() > 0) {

          int lab = partial_labels.back();
          partial_labels.pop_back();
          int rootx = find(lab, nodes);
          nodes[rootx] = find(ML, nodes);

        }

      }

      nextlabel = nextlabel + 1;

    }

  }

  // Rprintf("E = %u\n", nextlabel);

  // Second pass
  for (int i = 0; i < image.size(); i++) {

    if ( labels[i] > 0) {

      labels[i] = (int) find( labels[i], nodes);

    }

  }

  // SOLO FALTA REORDENAR LAS COMPONENTES CONEXAS

  // Rprintf("F\n");

  free(neighbours_offsets);
  free(neighs);
  free(label_neighs);
  free(nodes);
  free(expanded_idx);

  return labels;

}

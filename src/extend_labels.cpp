#include "BinaryHeapLabel.h"
// #include "pointers_neighborhoods.h"
#include <memory.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

#define INF 10000000

enum state {DEAD, ACTIVE, INACTIVE, BOUNDARY, SEED};

// [[Rcpp::export]]
IntegerVector extend_labels(IntegerVector pIn, IntegerVector maskImage) {

  IntegerVector dims = pIn.attr("dim");
  int ndims = dims.size();

  IntegerVector L(pIn.size());
  L.attr("dim") = dims;

  double* distances = (double*) malloc(pIn.size() * sizeof(double));

  memcpy(L.begin(), pIn.begin(), L.size() * sizeof(int));

  int N = pIn.size();

  CBinaryHeapLabelMinSort<int, double, int> heap(N);
  unsigned char* state = new unsigned char[N];

  int plane_size = dims[0] * dims[1];
  int row_size = dims[0];
  int neighbor_idx[6] = {-1, 1, -row_size, row_size, -plane_size, plane_size};

  int* coords = (int*) malloc(ndims * sizeof(int));

  // (1) Put seed into the heap and set their states to active, others set to inactive
  // NOTE: for efficiency purpose, boudaries of the matrix will not be processed.
  memset(state, BOUNDARY, N);
  for(int k = 1; k < dims[2] - 1; k++) {

    for(int j = 1;  j < dims[1] - 1; j++) {

      int idx_base = k * plane_size + j * row_size;

      for(int i = 1, idx = idx_base + 1; i < dims[0] - 1;  i++, idx++) {

        int mylabel = pIn[idx];

        if (mylabel > 0) {

          heap.Insert(idx, 0, mylabel);

          distances[idx] = 0;
          state[idx] = SEED;

        } else {

          state[idx] = INACTIVE;

          if (maskImage[idx] == 0) {

            L[idx] = 0;
            state[idx] = DEAD;

          }
        }
      }
    }
  }

  // (2) Extract a point from the heap
  int idx;
  double weight;
  int label;
  while( heap.Extract(&idx, &weight, &label) ) {

    // (3) if the state of the point is active ...
    if( state[idx] == DEAD ) continue;

    L[idx] = label;

    // (4) calculate the distance from its dead neighbors by solving ||\nabla u|| = 1;
    int x_orig = idx % dims[0];
    int y_orig = ((idx - x_orig) / dims[0]) % dims[1];
    int z_orig = ((idx - x_orig) / dims[0] - y_orig) / dims[1];

    // Get actual neighbours
    coords[0] = x_orig;
    coords[1] = y_orig;
    coords[2] = z_orig;

    double a, b, c, v1, v2, delta;
    int idx1, idx2;
    if( state[idx] != SEED ) {

      a = b = 0;  c = -1;

      // X:
      if (coords[0] > 0) {

        int n_idx = idx + neighbor_idx[0];
        v1 = (state[n_idx] == DEAD) ? distances[n_idx] : INF;

      }

      if (coords[0] < dims[0] - 1) {

        int n_idx = idx + neighbor_idx[1];
        v2 = (state[n_idx] == DEAD) ? distances[n_idx] : INF;

      }

      if( v1 > v2 ) v1 = v2;
      if( v1 != INF ) { a += 1; b += -2 * v1; c += v1 * v1; }

      // Y:
      if (coords[1] > 0) {

        int n_idx = idx + neighbor_idx[2];
        v1 = (state[n_idx] == DEAD) ? distances[n_idx] : INF;

      }

      if (coords[1] < dims[1] - 1) {

        int n_idx = idx + neighbor_idx[3];
        v2 = (state[n_idx] == DEAD) ? distances[n_idx] : INF;

      }

      if( v1 > v2 ) v1 = v2;
      if( v1 != INF ) { a += 1; b += -2 * v1; c += v1 * v1; }

      // Z:
      if (coords[2] > 0) {

        int n_idx = idx + neighbor_idx[4];
        v1 = (state[n_idx] == DEAD) ? distances[n_idx] : INF;

      }

      if (coords[2] < dims[2] - 1) {

        int n_idx = idx + neighbor_idx[5];
        v2 = (state[n_idx] == DEAD) ? distances[n_idx] : INF;


      }

      if( v1 > v2 ) v1 = v2;
      if( v1 != INF ) { a += 1; b += -2 * v1; c += v1 * v1; }

      // NOTE: a shouldn't be 0!
      delta = b * b - 4 * a * c;

      distances[idx] = ( delta > 0 ) ? .5 * (-b + sqrt(delta)) / a : 0;

    } else distances[idx] = 0;

    // (5) set the state to dead
    state[idx] = DEAD;

    // (6) search the 6-neighborhoods for non-deads (a) calculate the distance (b) push them to the heap
    // (c) set status to active
    // Add neighbours in three axis:
    // X:
    if (coords[0] > 0) {

      int n_idx = idx + neighbor_idx[0];

      if( state[n_idx] == INACTIVE || state[n_idx] == ACTIVE ) {

        heap.Insert(n_idx, distances[idx] + 1, label);
        state[n_idx] = ACTIVE;

      }

    }

    if (coords[0] < dims[0] - 1) {

      int n_idx = idx + neighbor_idx[1];

      if( state[n_idx] == INACTIVE || state[n_idx] == ACTIVE ) {

        heap.Insert(n_idx, distances[idx] + 1, label);
        state[n_idx] = ACTIVE;

      }

    }

    // Y:
    if (coords[1] > 0) {

      int n_idx = idx + neighbor_idx[2];

      if( state[n_idx] == INACTIVE || state[n_idx] == ACTIVE ) {

        heap.Insert(n_idx, distances[idx] + 1, label);
        state[n_idx] = ACTIVE;

      }

    }

    if (coords[1] < dims[1] - 1) {

      int n_idx = idx + neighbor_idx[3];

      if( state[n_idx] == INACTIVE || state[n_idx] == ACTIVE ) {

        heap.Insert(n_idx, distances[idx] + 1, label);
        state[n_idx] = ACTIVE;

      }

    }

    // Z:
    if (coords[2] > 0) {

      int n_idx = idx + neighbor_idx[4];

      if( state[n_idx] == INACTIVE || state[n_idx] == ACTIVE ) {

        heap.Insert(n_idx, distances[idx] + 1, label);
        state[n_idx] = ACTIVE;

      }

    }

    if (coords[2] < dims[2] - 1) {

      int n_idx = idx + neighbor_idx[5];

      if( state[n_idx] == INACTIVE || state[n_idx] == ACTIVE ) {

        heap.Insert(n_idx, distances[idx] + 1, label);
        state[n_idx] = ACTIVE;

      }

    }
  }

  delete [] state;
  free(coords);

  free(distances);

  return L;

}

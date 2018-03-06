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
NumericVector get_windows_at(NumericVector V, int width, IntegerVector x, IntegerVector y, IntegerVector z) {

  // if (width % 2 == 0) width++;
  int real_width = width;
  int disp = 0;
  if (width % 2 == 0) {

    real_width++;
    disp = 1;

  }

  int radius = (real_width + 1) / 2;

  IntegerVector dims = V.attr("dim");
  int n_volumes = 1;

  if (dims.size() >= 4) {

    n_volumes = dims[3];

  }

  int n_neighbours = pow(width, 3) * n_volumes;

  // Count how many windows we'll have
  int count = 0;

  NumericMatrix res(x.size(), n_neighbours + 3);

  // For every seed voxel, compute the window
  for (int i = 0; i < x.size(); i++) {

    res(count, 0) = x[i];
    res(count, 1) = y[i];
    res(count, 2) = z[i];

    int inner_count = 0;

    for (int volume = 0; volume < n_volumes; volume++) {

      for (int dz = -radius + 1; dz < radius - disp; dz++) {

        for (int dy = -radius + 1; dy < radius - disp; dy++) {

          for (int dx = -radius + 1; dx < radius - disp; dx++) {

            if ((x[i] + dx >= 0) & (x[i] + dx < dims[0]) & (y[i] + dy >= 0) & (y[i] + dy < dims[1]) & (z[i] + dz >= 0) & (z[i] + dz < dims[2])) {

              int offset = (x[i] + dx) + dims[0] * (y[i] + dy) + dims[0] * dims[1] * (z[i] + dz);
              offset += volume * dims[0] * dims[1] * dims[2];

              res(count, inner_count + 3) = V[offset];

            }

            inner_count++;

          }

        }

      }

    }

    count++;

  }

  return res;

}


// [[Rcpp::export]]
void results_to_volume(NumericVector V,
                       int width,
                       NumericVector res,
                       NumericVector counts,
                       IntegerVector x, IntegerVector y, IntegerVector z) {

  // if (width % 2 == 0) width++;
  int real_width = width;
  int disp = 0;
  if (width % 2 == 0) {

    real_width++;
    disp = 1;

  }

  int radius = (real_width + 1) / 2;

  IntegerVector dims = V.attr("dim");
  IntegerVector target_dims = res.attr("dim");

  for (int i = 0; i < x.size(); i++) {

    int inner_count = 0;

    for (int dz = -radius + 1; dz < radius - disp; dz++) {

      for (int dy = -radius + 1; dy < radius - disp; dy++) {

        for (int dx = -radius + 1; dx < radius - disp; dx++) {

          if ((x[i] + dx >= 0) & (x[i] + dx < target_dims[0]) & (y[i] + dy >= 0) & (y[i] + dy < target_dims[1]) & (z[i] + dz >= 0) & (z[i] + dz < target_dims[2])) {

            int offset = (x[i] + dx) + target_dims[0] * (y[i] + dy) + target_dims[0] * target_dims[1] * (z[i] + dz);

            res[offset] += V[i + dims[0] * inner_count];
            counts[offset] += 1;

          }

          inner_count++;

        }

      }

    }

  }

}

// [[Rcpp::export]]
void results_to_volume_label(NumericVector V,
                             int width,
                             NumericVector res,
                             IntegerVector x, IntegerVector y, IntegerVector z) {

  // if (width % 2 == 0) width++;
  int real_width = width;
  int disp = 0;
  if (width % 2 == 0) {

    real_width++;
    disp = 1;

  }

  int radius = (real_width + 1) / 2;

  IntegerVector dims = V.attr("dim");
  IntegerVector target_dims = res.attr("dim");

  // Rprintf("Vector dims = (%u,%u)\n", dims[0], dims[1]);

  for (int i = 0; i < x.size(); i++) {

    int inner_count = 0;

    for (int dz = -radius + 1; dz < radius - disp; dz++) {

      for (int dy = -radius + 1; dy < radius - disp; dy++) {

        for (int dx = -radius + 1; dx < radius - disp; dx++) {

          int offset = (x[i] + dx) + target_dims[0] * (y[i] + dy) + target_dims[0] * target_dims[1] * (z[i] + dz);

          if (V[i + dims[0] * inner_count] > res[offset])
            res[offset] = V[i + dims[0] * inner_count];
          inner_count++;

        }

      }

    }

  }

}

// [[Rcpp::export]]
void results_to_volume_label_with_distance(NumericVector V,
                                           int width,
                                           NumericVector res,
                                           NumericVector last_distance,
                                           IntegerVector x, IntegerVector y, IntegerVector z) {

  // if (width % 2 == 0) width++;
  int real_width = width;
  int disp = 0;
  if (width % 2 == 0) {

    real_width++;
    disp = 1;

  }

  int radius = (real_width + 1) / 2;

  IntegerVector dims = V.attr("dim");
  IntegerVector target_dims = res.attr("dim");

  // Rprintf("Vector dims = (%u,%u)\n", dims[0], dims[1]);

  for (int i = 0; i < x.size(); i++) {

    int inner_count = 0;

    for (int dz = -radius + 1; dz < radius - disp; dz++) {

      for (int dy = -radius + 1; dy < radius - disp; dy++) {

        for (int dx = -radius + 1; dx < radius - disp; dx++) {

          if ((x[i] + dx >= 0) & (x[i] + dx < target_dims[0]) & (y[i] + dy >= 0) & (y[i] + dy < target_dims[1]) & (z[i] + dz >= 0) & (z[i] + dz < target_dims[2])) {

            int offset = (x[i] + dx) + target_dims[0] * (y[i] + dy) + target_dims[0] * target_dims[1] * (z[i] + dz);

            if (dx * dx + dy * dy + dz * dz < last_distance[offset]) {

              res[offset] = V[i + dims[0] * inner_count];
              last_distance[offset] = dx * dx + dy * dy + dz * dz;

            }

          }

          inner_count++;

        }

      }

    }

  }

}

#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
#include "similarity_measures.h"
#include "utils.h"
#include "neighbours.h"

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

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

int count_eligible_voxels(int ndims,
                          int* dims,
                          int* mask,
                          int stride,
                          int* voxel_lookup_table) {

  int n_voxels_image = 1;
  for (int i = 0; i < ndims; i++) {

    n_voxels_image *= dims[i];

  }

  int count = 0;

  for (int x = 0; x < dims[0]; x++) {

    for (int y = 0; y < dims[1]; y++) {

      for (int z = 0; z < dims[2]; z++) {

        int i = x + dims[0] * y + dims[0] * dims[1] * z;

        if (mask[i] > 0 &&
            x >= 0 && x < dims[0] &&
            y >= 0 && y < dims[1] &&
            z >= 0 && z < dims[2] &&
            x % stride == 0 &&
            y % stride == 0 &&
            z % stride == 0) {

          voxel_lookup_table[i] = count;

          count++;

        } else {

          voxel_lookup_table[i] = -1;

        }


      }

    }

  }

  return count;

}


// [[Rcpp::export]]
int count_elegible(NumericVector image,
                   int patch_size,
                   int search_size,
                   int stride,
                   IntegerVector voxel_lookup_table) {

  IntegerVector dims = image.attr("dim");

  int* mask = (int*)malloc(image.size() * sizeof(int));
  for (int i = 0; i < image.size(); i++) {

    mask[i] = 1;

  }

  // Count elegible voxels
  int actual_voxels = count_eligible_voxels(3,
                                            dims.begin(),
                                            mask,
                                            stride,
                                            voxel_lookup_table.begin());

  free(mask);

  return actual_voxels;

}

// [[Rcpp::export]]
int count_elegible_masked(NumericVector image,
                          IntegerVector mask,
                          int patch_size,
                          int search_size,
                          int stride,
                          IntegerVector voxel_lookup_table) {

  IntegerVector dims = image.attr("dim");

  // Count elegible voxels
  int actual_voxels = count_eligible_voxels(3,
                                            dims.begin(),
                                            mask.begin(),
                                            stride,
                                            voxel_lookup_table.begin());

  return actual_voxels;

}


// [[Rcpp::export]]
void constrained_initialization_omp(NumericVector input_image,
                                    NumericVector template4D,
                                    int patch_size,
                                    int search_size,
                                    int actual_voxels,
                                    IntegerVector voxel_lookup_table,
                                    IntegerVector kANN,
                                    int ncores = 1) {


  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");
  int n_templates = template_dims[3];

  // Initialize similarities etc
  // First, indices for patch neighbours and search neighbours

  if (search_size % 2 == 0) search_size++;

  // int limits = (int)((patch_size + search_size) / 2 + 1);

  // int n_voxels_image = input_image.size();

  search_size = search_size / 2;
  if (search_size % 2 == 0) search_size++;

  // Added an omp pragma directive to parallelize the loop with ncores
  // #pragma omp parallel for num_threads(ncores)
  for (int x = 0; x < dims[0]; x++) {

    for (int y = 0; y < dims[1]; y++) {

      for (int z = 0; z < dims[2]; z++) {

        int i = x + dims[0] * y + dims[0] * dims[1] * z;
        int idx = voxel_lookup_table[i];

        if (idx >= 0) {

          for (int k1 = 0; k1 < n_templates; k1++) {

            int candidate;
            int x_displacement, y_displacement, z_displacement;

            do {

              x_displacement = rand_in_range(-search_size, search_size);

            } while ((x + x_displacement < 0) | (x + x_displacement >= dims[0]));

            do {

              y_displacement = rand_in_range(-search_size, search_size);

            } while ((y + y_displacement < 0) | (y + y_displacement >= dims[1]));

            do {

              z_displacement = rand_in_range(-search_size, search_size);

            } while ((z + z_displacement < 0) | (z + z_displacement >= dims[2] ));

            int absolute_displacement = x_displacement + dims[0] * y_displacement + dims[0] * dims[1] * z_displacement;
            candidate = i + absolute_displacement;

            kANN[k1 * actual_voxels + idx] = candidate;

          }

        }


      }

    }

  }

}

// [[Rcpp::export]]
void all_patches_similarity_omp(NumericVector input_image,
                                NumericVector template4D,
                                int actual_voxels,
                                IntegerVector voxel_lookup_table,
                                IntegerVector patch_neighbours,
                                IntegerVector kANN,
                                NumericVector similarities,
                                int method = 0,
                                int ncores = 1) {


  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");
  int n_templates = template_dims[3];

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(guided)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    int idx = voxel_lookup_table[voxel];

    if (idx >=0) {

      input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
      input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

      temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
      temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

      get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
      get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch, 0, n_voxels_image);
      if (method == 0) normalize(input_values, n_neighbours_patch);

      for (int k1 = 0; k1 < n_templates; k1++) {

        int candidate = k1 * n_voxels_image + kANN[k1 * actual_voxels + idx];

        get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

        get_image_value(template4D.begin(),
                        temp_neighs,
                        temp_values,
                        n_neighbours_patch,
                        k1 * n_voxels_image,
                        (k1 + 1) * n_voxels_image);

        double match = similarity(input_values, temp_values, n_neighbours_patch, method);

        similarities[k1 * actual_voxels + idx] = match;

      }

      free(input_values);
      free(temp_values);
      free(input_patch);
      free(temp_neighs);

    }


  }


}

// [[Rcpp::export]]
void propagation_step_omp(NumericVector input_image,
                          NumericVector template4D,
                          int actual_voxels,
                          IntegerVector voxel_lookup_table,
                          IntegerVector patch_neighbours,
                          IntegerVector kANN,
                          int direction,
                          int patch_size,
                          int stride,
                          NumericVector similarities,
                          int method = 0,
                          int ncores = 1) {

  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");
  int n_templates = template_dims[3];

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int* neighbour_offset = (int*) malloc(3 * sizeof(int));
  neighbour_offset[0] = direction * stride;
  neighbour_offset[1] = neighbour_offset[0] * dims[0];
  neighbour_offset[2] = neighbour_offset[1] * dims[1];

  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(guided)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    int idx = voxel_lookup_table[voxel];

    if (idx >= 0) {

      input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
      input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

      temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
      temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

      // Normalize input patch
      get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
      get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch, 0, n_voxels_image);
      if (method == 0) normalize(input_values, n_neighbours_patch);

      // Loop over all k ANN.
      for (int k1 = 0; k1 < n_templates; k1++) {

        // Loop over all 6 neighbours
        for (int neighbour = 0; neighbour < 3; neighbour++) {

          int loc = voxel + neighbour_offset[neighbour];

          if (loc < 0 || loc >= n_voxels_image) continue;

          int idx_loc = voxel_lookup_table[loc];

          if (idx_loc >= 0) {

            int original_candidate = kANN[idx_loc + k1 * actual_voxels];

            int candidate = original_candidate - neighbour_offset[neighbour];

            int relative_index = candidate;
            int x_orig = relative_index % dims[0];
            int y_orig = ((relative_index - x_orig) / dims[0]) % dims[1];
            int z_orig = ((relative_index - x_orig) / dims[0] - y_orig) / dims[1];

            if (x_orig < patch_size || x_orig > dims[0] - patch_size ||
                y_orig < patch_size || y_orig > dims[1] - patch_size ||
                z_orig < patch_size || z_orig > dims[2] - patch_size)
              continue;

            get_neighbours(k1 * n_voxels_image + candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

            get_image_value(template4D.begin(),
                            temp_neighs,
                            temp_values,
                            n_neighbours_patch,
                            k1 * n_voxels_image,
                            (k1 + 1) * n_voxels_image);

            double match = similarity(input_values, temp_values, n_neighbours_patch, method);

            if (match < similarities[k1 * actual_voxels + idx]) {

              similarities[k1 * actual_voxels + idx] = match;
              kANN[k1 * actual_voxels + idx] = candidate;

            }

          }

        }

      }

      free(input_values);
      free(temp_values);
      free(input_patch);
      free(temp_neighs);

    }



  }

  free(neighbour_offset);

}



// [[Rcpp::export]]
void constrained_random_search_omp(NumericVector input_image,
                                   NumericVector template4D,
                                   int actual_voxels,
                                   IntegerVector voxel_lookup_table,
                                   IntegerVector kANN,
                                   int patch_size,
                                   IntegerVector patch_neighbours,
                                   int search_size_max,
                                   NumericVector similarities,
                                   int max_random_neighbours,
                                   int method = 0,
                                   int ncores = 1) {


  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");

  IntegerVector template_dims = template4D.attr("dim");
  int n_templates = template_dims[3];

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int search_size = search_size_max / 2;
  int lower_limit = 0; //patch_size % 2 == 0 ? (patch_size + 3) / 2 : (patch_size + 1) / 2;

  while (search_size > 1)  {

    if (search_size % 2 == 0) search_size++;

    // Rprintf("search_size = %u\n", search_size);

    // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(guided)
    for (int voxel = 0; voxel < n_voxels_image; voxel++) {

      int idx_voxel = voxel_lookup_table[voxel];

      if (idx_voxel >= 0) {

        input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
        input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

        temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
        temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

        // Normalize input patch
        get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
        get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch, 0, n_voxels_image);
        if (method == 0) normalize(input_values, n_neighbours_patch);

        for (int k1 = 0; k1 < n_templates; k1++) {

          int idx = kANN[idx_voxel + k1 * actual_voxels];

          int relative_index = idx;

          int x_orig = relative_index % dims[0];
          int y_orig = ((relative_index - x_orig) / dims[0]) % dims[1];
          int z_orig = ((relative_index - x_orig) / dims[0] - y_orig) / dims[1];

          for (int neigh = 0; neigh < max_random_neighbours; neigh++) {

            int x_displacement, y_displacement, z_displacement;
            do {

              x_displacement = rand_in_range(-search_size, search_size);

            } while ((x_orig + x_displacement < lower_limit) || (x_orig + x_displacement >= dims[0] - lower_limit));

            do {

              y_displacement = rand_in_range(-search_size, search_size);

            } while ((y_orig + y_displacement < lower_limit) || (y_orig + y_displacement >= dims[1] - lower_limit));

            do {

              z_displacement = rand_in_range(-search_size, search_size);

            } while ((z_orig + z_displacement < lower_limit) || (z_orig + z_displacement >= dims[2] - lower_limit));

            int absolute_displacement = x_displacement + dims[0] * y_displacement + dims[0] * dims[1] * z_displacement;
            int candidate = idx + absolute_displacement;

            get_neighbours(k1 * n_voxels_image + candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

            get_image_value(template4D.begin(),
                            temp_neighs,
                            temp_values,
                            n_neighbours_patch,
                            k1 * n_voxels_image,
                            (k1 + 1) * n_voxels_image);

            double match = similarity(input_values, temp_values, n_neighbours_patch, method);

            if (match < similarities[k1 * actual_voxels + idx_voxel]) {

              similarities[k1 * actual_voxels + idx_voxel] = match;
              kANN[k1 * actual_voxels + idx_voxel] = candidate;

            }

          }

        }

        free(input_values);
        free(temp_values);
        free(input_patch);
        free(temp_neighs);

      }

    }

    search_size = search_size - 2;

  }

}


// // [[Rcpp::export]]
// void label_fusion_omp(IntegerVector labels4D,
//                       int actual_voxels,
//                       IntegerVector voxel_lookup_table,
//                       IntegerVector label_ids,
//                       IntegerVector kANN,
//                       IntegerVector patch_neighbours,
//                       int k,
//                       double lambda,
//                       double sigma2,
//                       NumericVector match,
//                       NumericVector new_voting,
//                       int ncores = 2) {
//
//   int n_neighbours_patch = patch_neighbours.size();
//   int n_voxels_image = voxel_lookup_table.size();
//   IntegerVector dims = labels4D.attr("dim");
//   int n_labels = label_ids.size();
//
//   double* similarities;
//   int* labels;
//   double* distances;
//
//   int* patch_seg;
//
//   int* input_patch;
//
//   int* temp_neighs;
//
//   int* counts = (int*) malloc(n_voxels_image * sizeof(int));
//   memset(counts, 0, n_voxels_image * sizeof(int));
//
//   int idx;
//
//   for (int voxel = 0; voxel < n_voxels_image; voxel++) {
//
//     idx = voxel_lookup_table[voxel];
//
//     if (idx >= 0) {
//
//       similarities = (double*) malloc(n_neighbours_patch * k * sizeof(double));
//       labels = (int*) malloc(n_neighbours_patch * k * sizeof(int));
//       distances = (double*) malloc(n_neighbours_patch * k * sizeof(double));
//
//       input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
//
//       // Compute the graylevel of the input image in the patch around given voxel
//       get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
//
//       // Added an omp pragma directive to parallelize the loop with ncores
// #pragma omp parallel for num_threads(ncores) private(temp_neighs, patch_seg) schedule(static)
//       for (int k1 = 0; k1 < k; k1++) {
//
//         int best_match = kANN[k1 * actual_voxels + idx];
//         int relative_index = best_match % n_voxels_image;
//
//         double distance = 0.0;
//
//         if (sigma2 > 0) {
//
//           int x_voxel = voxel % dims[0];
//           int x_candidate = relative_index % dims[0];
//           int y_voxel = ((voxel - x_voxel) / dims[0]) % dims[1];
//           int y_candidate = ((relative_index - x_candidate) / dims[0]) % dims[1];
//           int z_voxel = ((voxel - x_voxel) / dims[0] - y_voxel) / dims[1];
//           int z_candidate = ((relative_index - x_candidate) / dims[0] - y_candidate) / dims[1];
//
//           distance = (x_voxel - x_candidate) * (x_voxel - x_candidate);
//           distance += (y_voxel - y_candidate) * (y_voxel - y_candidate);
//           distance += (z_voxel - z_candidate) * (z_voxel - z_candidate);
//
//         }
//
//         temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
//         patch_seg = (int*) malloc(n_neighbours_patch * sizeof(int));
//
//         get_neighbours(kANN[k1 * actual_voxels + idx], patch_neighbours.begin(), temp_neighs, n_neighbours_patch);
//         get_label_value(labels4D.begin(), temp_neighs, patch_seg, n_neighbours_patch);
//
//         double simil = match[k1 * actual_voxels + idx];
//
//         // If we use the distances to neighboring patches in the search to weight the similarity
//         for (int neigh = 0; neigh < n_neighbours_patch; neigh++) {
//
//           similarities[k1 * n_neighbours_patch + neigh] = simil;
//           labels[k1 * n_neighbours_patch + neigh] = patch_seg[neigh];
//           distances[k1 * n_neighbours_patch + neigh] = distance;
//
//         }
//
//         free(temp_neighs);
//         free(patch_seg);
//
//       } // For all k
//
//
//       // Normalize the similarities
//       // First, compute the minimum similarity
//       double min_simil = 1.e100;
//       for (int temp = 0; temp < k; temp++) {
//
//         if (similarities[temp * n_neighbours_patch] < min_simil) {
//
//           min_simil = similarities[temp * n_neighbours_patch];
//
//         }
//
//       }
//       // double min_simil = find_minimum(similarities, k * n_neighbours_patch);
//
//
//       double h2 = lambda * min_simil + 1.e-15;
//       for (int i = 0; i < k * n_neighbours_patch; i++) {
//
//         similarities[i] = exp(-similarities[i] / h2 -  distances[i] / (sigma2 + 1.e-15));
//
//       }
//
//       // Compute the voting for each neighbour
//       double cumulative_similarity = 0;
//       for (int temp = 0; temp < k; temp++) {
//
//         cumulative_similarity += similarities[temp * n_neighbours_patch];
//
//       }
//
//       // For all labels except the first (probably label == 0)
// #pragma omp parallel for num_threads(ncores) schedule(dynamic)
//       for (int lb = 1; lb < n_labels; lb++) {
//
//         for (int i = 0; i < n_neighbours_patch; i++) {
//
//           int which_voxel = input_patch[i];
//
//           if (which_voxel > -1) {
//
//             counts[which_voxel]++;
//
//             double similarity = 0;
//
//             for (int temp = 0; temp < k; temp++) {
//
//               if (labels[temp * n_neighbours_patch + i] == label_ids[lb]) {
//
//                 similarity += similarities[temp * n_neighbours_patch + i];
//
//               }
//
//             }
//
//             new_voting[lb * n_voxels_image + which_voxel] += similarity / cumulative_similarity;
//
//           }
//
//         }
//
//       }
//
//       free(input_patch);
//       free(similarities);
//       free(labels);
//       free(distances);
//
//     }
//
//   }
//
// #pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
//   for (int voxel = 0; voxel < n_voxels_image; voxel++) {
//
//     double simil_0 = 1;
//     for (int lb = 1; lb < n_labels; lb++) {
//
//       simil_0 -= new_voting[lb * n_voxels_image + voxel];
//
//     }
//
//     new_voting[voxel] = simil_0;
//
//   }
//
// #pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
//   for (int voxel = 0; voxel < n_voxels_image; voxel++) {
//
//     if (counts[voxel] > 0) {
//
//       for (int lb = 0; lb < n_labels; lb++) {
//
//         new_voting[lb * n_voxels_image + voxel] /= counts[voxel];
//
//       }
//
//     }
//
//   }
//
//   free(counts);
//
// }


// [[Rcpp::export]]
void label_fusion2_omp(IntegerVector labels4D,
                       int actual_voxels,
                       IntegerVector voxel_lookup_table,
                       IntegerVector label_ids,
                       IntegerVector kANN,
                       IntegerVector patch_neighbours,
                       double lambda,
                       double sigma2,
                       NumericVector match,
                       NumericVector new_voting,
                       int ncores = 2) {

  int n_neighbours_patch = patch_neighbours.size();
  int n_voxels_image = voxel_lookup_table.size();
  IntegerVector dims = labels4D.attr("dim");
  int n_labels = label_ids.size();

  int n_templates = dims[3];

  int counts;
  int max_label = max(labels4D);

  int* new_labels = (int*) malloc((max_label + 1) * sizeof(int));

  for (int lb = 0; lb <= max_label; lb++) {

    new_labels[lb] = -1;

  }

  for (int lb = 0; lb < n_labels; lb++) {

    new_labels[label_ids[lb]] = lb;

  }

#pragma omp parallel for num_threads(ncores) schedule(guided)
  for (int idx = 0; idx < actual_voxels; idx++) {

    double min_simil = 1.e10;
    for (int k1 = 0; k1 < n_templates; k1++) {

      if (min_simil > match[k1 * actual_voxels + idx])
        min_simil = match[k1 * actual_voxels + idx];

    }

    double h2 = min_simil * lambda + 1.e-15;
    for (int k1 = 0; k1 < n_templates; k1++) {

      match[k1 * actual_voxels + idx] = exp( match[k1 * actual_voxels + idx] / h2 );

    }

    double cumulative = 0;
    for (int k1 = 0; k1 < n_templates; k1++) {

      cumulative += match[k1 * actual_voxels + idx];

    }

    for (int k1 = 0; k1 < n_templates; k1++) {

      match[k1 * actual_voxels + idx] /= cumulative;

    }

  }


#pragma omp parallel for num_threads(ncores) private(counts) schedule(guided)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    counts = 0;

    // Loop over all neighbours
    for (int neighbour = 0; neighbour < n_neighbours_patch; neighbour++) {

      int my_neighbour = voxel + patch_neighbours[neighbour];

      if ((my_neighbour < 0) || (my_neighbour >= n_voxels_image)) continue;

      // If idx_patch_center >= 0 is because it is the center of a patch
      int idx_patch_center = voxel_lookup_table[my_neighbour];

      if (idx_patch_center >= 0) {

        int relative_difference = voxel - my_neighbour;

        for (int k1 = 0; k1 < n_templates; k1++) {

          int candidate =  k1 * n_voxels_image + kANN[k1 * actual_voxels + idx_patch_center];
          int voxel_candidate = candidate + relative_difference;

          if (voxel_candidate < k1 * n_voxels_image) {

            voxel_candidate = k1 * n_voxels_image;

          }

          if (voxel_candidate >= (k1 + 1) * n_voxels_image) {

            voxel_candidate = (k1 + 1) * n_voxels_image - 1;
          }

          int lb = new_labels[labels4D[voxel_candidate]];

          if (lb >= 0)
            new_voting[lb * n_voxels_image + voxel] += match[k1 * actual_voxels + idx_patch_center];

        }


        counts++;

      } // Patch

    } // All neighbours


    if (counts > 0) {

      for (int lb = 1; lb < n_labels; lb++) {

        new_voting[lb * n_voxels_image + voxel] /= counts;

      }

    }

    double simil_0 = 1;
    for (int lb = 1; lb < n_labels; lb++) {

      simil_0 -= new_voting[lb * n_voxels_image + voxel];

    }

    new_voting[voxel] = simil_0;

  } // All voxels


  free(new_labels);

}


// [[Rcpp::export]]
void label_fusion_omp_fast(IntegerVector labels4D,
                           int actual_voxels,
                           IntegerVector voxel_lookup_table,
                           IntegerVector label_ids,
                           IntegerVector kANN,
                           IntegerVector patch_neighbours,
                           double lambda,
                           double sigma2,
                           NumericVector match,
                           IntegerVector new_voting,
                           NumericVector new_sim,
                           int ncores = 2) {

  int n_neighbours_patch = patch_neighbours.size();
  int n_voxels_image = voxel_lookup_table.size();
  IntegerVector dims = labels4D.attr("dim");
  int n_labels = label_ids.size();

  int n_templates = dims[3];

  int counts;
  int max_label = max(labels4D);

  int* new_labels = (int*) malloc((max_label + 1) * sizeof(int));

  for (int lb = 0; lb <= max_label; lb++) {

    new_labels[lb] = -1;

  }

  for (int lb = 0; lb < n_labels; lb++) {

    new_labels[label_ids[lb]] = lb;

  }

#pragma omp parallel for num_threads(ncores) schedule(guided)
  for (int idx = 0; idx < actual_voxels; idx++) {

    double min_simil = 1.e10;
    for (int k1 = 0; k1 < n_templates; k1++) {

      if (min_simil > match[k1 * actual_voxels + idx])
        min_simil = match[k1 * actual_voxels + idx];

    }

    double h2 = min_simil * lambda + 1.e-15;
    for (int k1 = 0; k1 < n_templates; k1++) {

      match[k1 * actual_voxels + idx] = exp( match[k1 * actual_voxels + idx] / h2 );

    }

    double cumulative = 0;
    for (int k1 = 0; k1 < n_templates; k1++) {

      cumulative += match[k1 * actual_voxels + idx];

    }

    for (int k1 = 0; k1 < n_templates; k1++) {

      match[k1 * actual_voxels + idx] /= cumulative;

    }

  }

#pragma omp parallel for num_threads(ncores) private(counts) schedule(guided)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    counts = 0;

    // Loop over all neighbours
    for (int neighbour = 0; neighbour < n_neighbours_patch; neighbour++) {

      int my_neighbour = voxel + patch_neighbours[neighbour];

      if ((my_neighbour < 0) || (my_neighbour >= n_voxels_image)) continue;

      // If idx_patch_center >= 0 is because it is the center of a patch
      int idx_patch_center = voxel_lookup_table[my_neighbour];

      if (idx_patch_center >= 0) {

        int relative_difference = voxel - my_neighbour;

        for (int k1 = 0; k1 < n_templates; k1++) {

          int candidate =  k1 * n_voxels_image + kANN[k1 * actual_voxels + idx_patch_center];
          int voxel_candidate = candidate + relative_difference;

          if (voxel_candidate < k1 * n_voxels_image) {

            voxel_candidate = k1 * n_voxels_image;

          }

          if (voxel_candidate >= (k1 + 1) * n_voxels_image) {

            voxel_candidate = (k1 + 1) * n_voxels_image - 1;
          }

          int lb = new_labels[labels4D[voxel_candidate]];

          if (lb >= 0)
            if (new_sim[k1 * n_voxels_image + voxel] < match[k1 * actual_voxels + idx_patch_center]) {

              new_voting[k1 * n_voxels_image + voxel] = lb;
              new_sim[k1 * n_voxels_image + voxel] = match[k1 * actual_voxels + idx_patch_center];

            }

        }

      } // Patch

    } // All neighbours

  } // All voxels

  free(new_labels);

}

// [[Rcpp::export]]
void label_fusion_mode(IntegerVector my_labels,
                       IntegerVector result,
                       NumericVector my_similarities,
                       int ncores = 1) {

  int n_voxels_image = result.size();

  Rprintf("n_voxels_image = %u\n", n_voxels_image);

  IntegerVector dims = my_labels.attr("dim");

  int n_templates = dims[3];

  int max_label = max(my_labels);

  double* new_labels = (double*) malloc((max_label + 1) * sizeof(double));

  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    // Rprintf("voxel = %u\n", voxel);

    for (int i = 0; i < max_label + 1; i++) {

      new_labels[i] = 0;

    }

    // Rprintf("  Initialized\n");

    for (int k = 0; k < n_templates; k++) {

      int lb = my_labels[k * n_voxels_image + voxel];

      new_labels[lb] += my_similarities[k * n_voxels_image + voxel];

    }

    // Rprintf("  Tabulated\n");

    double M = -1;
    int M_label;
    for (int i = 0; i < max_label + 1; i++) {

      if (M < new_labels[i]) {

        M = new_labels[i];
        M_label = i;

      }

    }

    // Rprintf("  Max\n");

    result[voxel] = M_label;

    // Rprintf("  Assigned\n");

  }

  // Rprintf("Out\n");

  free(new_labels);

}

// [[Rcpp::export]]
void image_fusion_omp(NumericVector labels4D,
                      int actual_voxels,
                      IntegerVector voxel_lookup_table,
                      IntegerVector kANN,
                      IntegerVector patch_neighbours,
                      double lambda,
                      double sigma2,
                      NumericVector match,
                      NumericVector new_voting,
                      IntegerVector voxel_candidate,
                      int ncores = 2) {

  int n_neighbours_patch = patch_neighbours.size();
  int n_voxels_image = voxel_lookup_table.size();
  IntegerVector dims = labels4D.attr("dim");

  int n_templates = dims[3];

  double counts;


#pragma omp parallel for num_threads(ncores)
  for (int idx = 0; idx < actual_voxels; idx++) {

    double min_simil = 1.e10;
    for (int k1 = 0; k1 < n_templates; k1++) {

      if (min_simil > match[k1 * actual_voxels + idx])
        min_simil = match[k1 * actual_voxels + idx];

    }

    double h2 = min_simil * lambda + 1.e-15;
    for (int k1 = 0; k1 < n_templates; k1++) {

      match[k1 * actual_voxels + idx] = exp( match[k1 * actual_voxels + idx] / h2 );

    }

    double cumulative = 0;
    for (int k1 = 0; k1 < n_templates; k1++) {

      cumulative += match[k1 * actual_voxels + idx];

    }

    for (int k1 = 0; k1 < n_templates; k1++) {

      match[k1 * actual_voxels + idx] /= cumulative;

    }

  }


#pragma omp parallel for num_threads(ncores) private(counts) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    counts = 0;
    voxel_candidate[voxel] = voxel;

    // Loop over all neighbours
    for (int neighbour = 0; neighbour < n_neighbours_patch; neighbour++) {

      int my_neighbour = voxel + patch_neighbours[neighbour];

      if ((my_neighbour < 0) || (my_neighbour >= n_voxels_image)) continue;

      // If idx_patch_center >= 0 is because it is the center of a patch
      int idx_patch_center = voxel_lookup_table[my_neighbour];

      if (idx_patch_center >= 0) {

        int relative_difference = voxel - my_neighbour;

        for (int k1 = 0; k1 < n_templates; k1++) {

          int candidate =  k1 * n_voxels_image + kANN[k1 * actual_voxels + idx_patch_center];
          voxel_candidate[voxel] = candidate + relative_difference;

          if (voxel_candidate[voxel] < k1 * n_voxels_image) {

            voxel_candidate[voxel] = k1 * n_voxels_image;

          }

          if (voxel_candidate[voxel] >= (k1 + 1) * n_voxels_image) {

            voxel_candidate[voxel] = (k1 + 1) * n_voxels_image - 1;

          }

          double lb = labels4D[voxel_candidate[voxel]];

          new_voting[voxel] += lb * match[k1 * actual_voxels + idx_patch_center];


          counts += match[k1 * actual_voxels + idx_patch_center];

        }

        // counts++;

      } // Patch

    } // All neighbours


    if (counts > 0) {

      new_voting[voxel] /= counts;

    }


  } // All voxels


}


// // [[Rcpp::export]]
// void label_fusion3_omp(IntegerVector labels4D,
//                        int actual_voxels,
//                        IntegerVector voxel_lookup_table,
//                        IntegerVector label_ids,
//                        IntegerVector kANN,
//                        IntegerVector patch_neighbours,
//                        int k,
//                        double lambda,
//                        double sigma2,
//                        NumericVector match,
//                        NumericVector new_voting,
//                        int ncores = 2) {
//
//   int n_neighbours_patch = patch_neighbours.size();
//   int n_voxels_image = voxel_lookup_table.size();
//   IntegerVector dims = labels4D.attr("dim");
//   int n_labels = label_ids.size();
//
//   double* similarities;
//   int* labels;
//   double* distances;
//
//   int* patch_seg;
//
//   int* input_patch;
//
//   int* temp_neighs;
//   int patch_size = (int)pow(n_neighbours_patch, 1/3);
//
//   int* counts = (int*) malloc(n_voxels_image * sizeof(int));
//   memset(counts, 0, n_voxels_image * sizeof(int));
//
//   int idx;
//
//   // Added an omp pragma directive to parallelize the loop with ncores
// #pragma omp parallel for num_threads(ncores) private(temp_neighs, patch_seg, similarities, labels, distances, input_patch) schedule(dynamic)
//   for (int init = 0; init < patch_size; init++) {
//
//     for (int x = init; x < dims[0]; x+=patch_size) {
//
//       for (int y = init; y < dims[1]; y+=patch_size) {
//
//         for (int z = init; z < dims[2]; z+=patch_size) {
//
//           int voxel = x + dims[0] * y + dims[0] * dims[1] * z;
//
//           idx = voxel_lookup_table[voxel];
//
//           if (idx >= 0) {
//
//             similarities = (double*) malloc(n_neighbours_patch * k * sizeof(double));
//             labels = (int*) malloc(n_neighbours_patch * k * sizeof(int));
//             distances = (double*) malloc(n_neighbours_patch * k * sizeof(double));
//
//             input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
//
//             // Compute the graylevel of the input image in the patch around given voxel
//             get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
//
//             for (int k1 = 0; k1 < k; k1++) {
//
//               int best_match = kANN[k1 * actual_voxels + idx];
//               int relative_index = best_match % n_voxels_image;
//
//               double distance = 0.0;
//
//               if (sigma2 > 0) {
//
//                 int x_voxel = voxel % dims[0];
//                 int x_candidate = relative_index % dims[0];
//                 int y_voxel = ((voxel - x_voxel) / dims[0]) % dims[1];
//                 int y_candidate = ((relative_index - x_candidate) / dims[0]) % dims[1];
//                 int z_voxel = ((voxel - x_voxel) / dims[0] - y_voxel) / dims[1];
//                 int z_candidate = ((relative_index - x_candidate) / dims[0] - y_candidate) / dims[1];
//
//                 distance = (x_voxel - x_candidate) * (x_voxel - x_candidate);
//                 distance += (y_voxel - y_candidate) * (y_voxel - y_candidate);
//                 distance += (z_voxel - z_candidate) * (z_voxel - z_candidate);
//
//               }
//
//               temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
//               patch_seg = (int*) malloc(n_neighbours_patch * sizeof(int));
//
//               get_neighbours(kANN[k1 * actual_voxels + idx], patch_neighbours.begin(), temp_neighs, n_neighbours_patch);
//               get_label_value(labels4D.begin(), temp_neighs, patch_seg, n_neighbours_patch);
//
//               double simil = match[k1 * actual_voxels + idx];
//
//               // If we use the distances to neighboring patches in the search to weight the similarity
//               for (int neigh = 0; neigh < n_neighbours_patch; neigh++) {
//
//                 similarities[k1 * n_neighbours_patch + neigh] = simil;
//                 labels[k1 * n_neighbours_patch + neigh] = patch_seg[neigh];
//                 distances[k1 * n_neighbours_patch + neigh] = distance;
//
//               }
//
//               free(temp_neighs);
//               free(patch_seg);
//
//             } // For all k
//
//
//             // Normalize the similarities
//             // First, compute the minimum similarity
//             double min_simil = 1.e100;
//             for (int temp = 0; temp < k; temp++) {
//
//               if (similarities[temp * n_neighbours_patch] < min_simil) {
//
//                 min_simil = similarities[temp * n_neighbours_patch];
//
//               }
//
//             }
//             // double min_simil = find_minimum(similarities, k * n_neighbours_patch);
//
//
//             double h2 = lambda * min_simil + 1.e-15;
//             for (int i = 0; i < k * n_neighbours_patch; i++) {
//
//               similarities[i] = exp(-similarities[i] / h2 -  distances[i] / (sigma2 + 1.e-15));
//
//             }
//
//             // Compute the voting for each neighbour
//             double cumulative_similarity = 0;
//             for (int temp = 0; temp < k; temp++) {
//
//               cumulative_similarity += similarities[temp * n_neighbours_patch];
//
//             }
//
//             // For all labels except the first (probably label == 0)
//             // #pragma omp parallel for num_threads(ncores) schedule(dynamic)
//             for (int lb = 1; lb < n_labels; lb++) {
//
//               for (int i = 0; i < n_neighbours_patch; i++) {
//
//                 int which_voxel = input_patch[i];
//
//                 if (which_voxel > -1) {
//
//                   counts[which_voxel]++;
//
//                   double similarity = 0;
//
//                   for (int temp = 0; temp < k; temp++) {
//
//                     if (labels[temp * n_neighbours_patch + i] == label_ids[lb]) {
//
//                       similarity += similarities[temp * n_neighbours_patch + i];
//
//                     }
//
//                   }
//
//                   new_voting[lb * n_voxels_image + which_voxel] += similarity / cumulative_similarity;
//
//                 }
//
//               }
//
//             }
//
//             free(input_patch);
//             free(similarities);
//             free(labels);
//             free(distances);
//
//           }
//
//
//         }
//
//       }
//
//     }
//
//   }
//
// #pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
//   for (int voxel = 0; voxel < n_voxels_image; voxel++) {
//
//     double simil_0 = 1;
//     for (int lb = 1; lb < n_labels; lb++) {
//
//       simil_0 -= new_voting[lb * n_voxels_image + voxel];
//
//     }
//
//     new_voting[voxel] = simil_0;
//
//   }
//
// #pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
//   for (int voxel = 0; voxel < n_voxels_image; voxel++) {
//
//     if (counts[voxel] > 0) {
//
//       for (int lb = 0; lb < n_labels; lb++) {
//
//         new_voting[lb * n_voxels_image + voxel] /= counts[voxel];
//
//       }
//
//     }
//
//   }
//
//   free(counts);
//
// }

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
                          int limits,
                          int stride,
                          int* voxel_lookup_table) {

  // Rprintf("A\n");
  int n_voxels_image = 1;
  for (int i = 0; i < ndims; i++) {

    n_voxels_image *= dims[i];

  }
  // Rprintf("B\n");

  int count = 0;

  // Rprintf("C\n");

  for (int x = 0; x < dims[0]; x++) {

    for (int y = 0; y < dims[1]; y++) {

      for (int z = 0; z < dims[2]; z++) {

        int i = x + dims[0] * y + dims[0] * dims[1] * z;

        if (x >= limits && x < dims[0] - limits &&
            y >= limits && y < dims[1] - limits &&
            z >= limits && z < dims[2] - limits &&
            x % stride == 0 &&
            y % stride == 0 &&
            z % stride == 0) {

          // Rprintf("In = %u = %u\n", i, count);
          voxel_lookup_table[i] = count;

          count++;

        } else {
          // Rprintf("Out = %u\n", i);

          voxel_lookup_table[i] = -1;

        }


      }

    }

  }

  // count--;

  // Rprintf("D\n");

  return count;

}


// [[Rcpp::export]]
int count_elegible(NumericVector image,
                   int patch_size,
                   int search_size,
                   int stride,
                   IntegerVector voxel_lookup_table) {

  IntegerVector dims = image.attr("dim");

  // Rprintf("n_voxels_image = %u\n", n_voxels_image);

  // Rprintf("Lookup\n");

  int limits = (int)((patch_size + search_size) / 2 + 1);

  // Count elegible voxels
  int actual_voxels = count_eligible_voxels(3,
                                            dims.begin(),
                                            limits,
                                            stride,
                                            voxel_lookup_table.begin());

  // Rprintf("Init = %u\n", actual_voxels);

  // free(voxel_lookup_table);

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
                                    int k,
                                    int ncores = 1) {


  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");
  int n_templates = template_dims[3];

  // Initialize similarities etc
  // First, indices for patch neighbours and search neighbours
  if (patch_size % 2 == 0) patch_size++;
  int n_neighbours_patch = pow(patch_size, 3);
  int* patch_neighbours = (int*) malloc(n_neighbours_patch * sizeof(int));
  get_neighbours_indices(dims.begin(), 3, patch_size, patch_neighbours);

  if (search_size % 2 == 0) search_size++;

  int limits = (int)((patch_size + search_size) / 2 + 1);

  int n_voxels_image = input_image.size();

  search_size = search_size / 2;
  if (search_size % 2 == 0) search_size++;

  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores)
  for (int x = 0; x < dims[0]; x++) {

    for (int y = 0; y < dims[1]; y++) {

      for (int z = 0; z < dims[2]; z++) {

        int i = x + dims[0] * y + dims[0] * dims[1] * z;
        int idx = voxel_lookup_table[i];

        if (idx >= 0) {

          for (int k1 = 0; k1 < k; k1++) {

            int candidate;
            int template_id = rand_in_range(1, n_templates) - 1;
            int x_displacement, y_displacement, z_displacement;

            do {

              x_displacement = rand_in_range(-search_size, search_size);

            } while ((x + x_displacement <= limits) | (x + x_displacement > dims[0] - limits));

            do {

              y_displacement = rand_in_range(-search_size, search_size);

            } while ((y + y_displacement <= limits) | (y + y_displacement > dims[1] - limits));

            do {

              z_displacement = rand_in_range(-search_size, search_size);

            } while ((z + z_displacement <= limits) | (z + z_displacement > dims[2] - limits));

            int absolute_displacement = x_displacement + dims[0] * y_displacement + dims[0] * dims[1] * z_displacement;
            candidate = i + absolute_displacement;

            kANN[k1 * actual_voxels + idx] = template_id * n_voxels_image + candidate;

          }

        }


      }

    }

  }

}

// [[Rcpp::export]]
void all_patches_similarity_omp(NumericVector input_image,
                                NumericVector template4D,
                                int k,
                                int actual_voxels,
                                IntegerVector voxel_lookup_table,
                                IntegerVector patch_neighbours,
                                IntegerVector kANN,
                                NumericVector similarities,
                                int ncores = 1) {


  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  // int CHUNK_SIZE = n_voxels_image / ncores;

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
      get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch);
      normalize(input_values, n_neighbours_patch);

      for (int k1 = 0; k1 < k; k1++) {

        int candidate = kANN[k1 * actual_voxels + idx];
        // Rprintf("k1 = %u, n_voxels_image = %u, voxel = %u, Candidate = %u\n", k1, n_voxels_image, voxel, candidate);

        get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

        get_image_value(template4D.begin(), temp_neighs, temp_values, n_neighbours_patch);

        double matchSum1 = 0, matchSSQ1 = 0;

        double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);

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
void all_patches_similarity_omp2(NumericVector input_image,
                                 NumericVector template4D,
                                 int k,
                                 int actual_voxels,
                                 IntegerVector voxel_lookup_table,
                                 IntegerVector voxel_array_index,
                                 IntegerVector patch_neighbours,
                                 IntegerVector kANN,
                                 NumericVector similarities,
                                 int ncores = 1) {


  // int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");
  // int n_templates = template_dims[3];

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int CHUNK_SIZE = actual_voxels / ncores;

  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(dynamic, CHUNK_SIZE)
  for (int voxel_index = 0; voxel_index < actual_voxels; voxel_index++) {

    int voxel = voxel_array_index[voxel_index];
    int idx = voxel_lookup_table[voxel_index];

    input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
    input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

    temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
    temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

    get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
    get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch);
    normalize(input_values, n_neighbours_patch);

    for (int k1 = 0; k1 < k; k1++) {

      int candidate = kANN[k1 * actual_voxels + idx];

      get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

      get_image_value(template4D.begin(), temp_neighs, temp_values, n_neighbours_patch);

      double matchSum1 = 0, matchSSQ1 = 0;

      double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);

      similarities[k1 * actual_voxels + idx] = match;

    }

    free(input_values);
    free(temp_values);
    free(input_patch);
    free(temp_neighs);

  }


}


// [[Rcpp::export]]
void propagation_step_omp(NumericVector input_image,
                          NumericVector template4D,
                          int actual_voxels,
                          IntegerVector voxel_lookup_table,
                          IntegerVector patch_neighbours,
                          IntegerVector kANN,
                          int k,
                          int direction,
                          int patch_size,
                          int stride,
                          NumericVector similarities,
                          int ncores = 1) {

  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  IntegerVector template_dims = template4D.attr("dim");
  // int n_templates = template_dims[3];

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int* neighbour_offset = (int*) malloc(3 * sizeof(int));
  neighbour_offset[0] = direction * stride;
  neighbour_offset[1] = neighbour_offset[0] * dims[0];
  neighbour_offset[2] = neighbour_offset[1] * dims[1];

  // int changes = 0;

  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    int idx = voxel_lookup_table[voxel];

    if (idx >= 0) {

      input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
      input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

      temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
      temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

      // Normalize input patch
      get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
      get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch);
      normalize(input_values, n_neighbours_patch);

      // Loop over all k ANN.
      for (int k1 = 0; k1 < k; k1++) {

        // Loop over all 6 neighbours
        for (int neighbour = 0; neighbour < 3; neighbour++) {

          int loc = voxel + neighbour_offset[neighbour];

          if (loc < 0 || loc >= n_voxels_image) continue;

          int idx_loc = voxel_lookup_table[loc];

          if (idx_loc >= 0) {

            int original_candidate = kANN[idx_loc + k1 * actual_voxels];
            int original_template = original_candidate / n_voxels_image;

            int candidate = original_candidate - neighbour_offset[neighbour];
            int new_template = candidate / n_voxels_image;

            if (original_template != new_template) {

              continue;

            }

            int relative_index = candidate - new_template * n_voxels_image;
            int x_orig = relative_index % dims[0];
            int y_orig = ((relative_index - x_orig) / dims[0]) % dims[1];
            int z_orig = ((relative_index - x_orig) / dims[0] - y_orig) / dims[1];

            if (x_orig < patch_size || x_orig > dims[0] - patch_size ||
                y_orig < patch_size || y_orig > dims[1] - patch_size ||
                z_orig < patch_size || z_orig > dims[2] - patch_size)
              continue;

            get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

            get_image_value(template4D.begin(), temp_neighs, temp_values, n_neighbours_patch);

            double matchSum1 = 0, matchSSQ1 = 0;

            double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);

            // double match = SSDPatch(input_values, temp_values, n_neighbours_patch, previous);

            if (match < similarities[k1 * actual_voxels + idx]) {

              // // Rprintf("Cambio\n");
              // changes++;

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

  // free(result);
  free(neighbour_offset);

  // Rprintf("Changes: %u\n", changes);

}


// [[Rcpp::export]]
void constrained_random_search_omp(NumericVector input_image,
                                   NumericVector template4D,
                                   int actual_voxels,
                                   IntegerVector voxel_lookup_table,
                                   IntegerVector kANN,
                                   int k,
                                   int patch_size,
                                   IntegerVector patch_neighbours,
                                   int search_size_max,
                                   NumericVector similarities,
                                   int max_random_neighbours,
                                   int ncores = 1) {


  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int search_size = search_size_max / 2;

  while (search_size > 1)  {

    // Rprintf("search_size = %u\n", search_size);

    if (search_size % 2 == 0) search_size++;

    // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(dynamic, 10000)
    for (int voxel = 0; voxel < n_voxels_image; voxel++) {

      int idx_voxel = voxel_lookup_table[voxel];

      if (idx_voxel >= 0) {

        input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
        input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

        temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
        temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

        // Rprintf("voxel = %u, idx_voxel = %d\n", voxel, idx_voxel);

        // Normalize input patch
        get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
        get_image_value(input_image.begin(), input_patch, input_values, n_neighbours_patch);
        normalize(input_values, n_neighbours_patch);

        for (int k1 = 0; k1 < k; k1++) {

          int idx = kANN[idx_voxel + k1 * actual_voxels];

          int relative_index = idx % n_voxels_image;

          // Rprintf("k1 = %d, actual_voxels = %d, n_voxels_image = %d, idx = %d, relative_index = %d\n",
          //         k1, actual_voxels, n_voxels_image, idx, relative_index);

          int x_orig = relative_index % dims[0];
          int y_orig = ((relative_index - x_orig) / dims[0]) % dims[1];
          int z_orig = ((relative_index - x_orig) / dims[0] - y_orig) / dims[1];

          // Rprintf("Candidate x: %u, y: %u, z: %u\n", x_orig, y_orig, z_orig);

          int x_displacement, y_displacement, z_displacement;
          do {

            x_displacement = rand_in_range(-search_size, search_size);

          } while ((x_orig + x_displacement < patch_size) || (x_orig + x_displacement >= dims[0] - patch_size));

          do {

            y_displacement = rand_in_range(-search_size, search_size);

          } while ((y_orig + y_displacement < patch_size) || (y_orig + y_displacement >= dims[1] - patch_size));

          do {

            z_displacement = rand_in_range(-search_size, search_size);

          } while ((z_orig + z_displacement < patch_size) || (z_orig + z_displacement >= dims[2] - patch_size));

          // Rprintf("New Candidate x: %u, y: %u, z: %u\n", x_orig + x_displacement, y_orig + y_displacement, z_orig + z_displacement);

          int absolute_displacement = x_displacement + dims[0] * y_displacement + dims[0] * dims[1] * z_displacement;
          int candidate = idx + absolute_displacement;

          // Rprintf("voxel / idx / candidate / template / total = %u / %u / %u / %u / %u\n",
          //         voxel, idx, candidate % n_voxels_image, candidate % n_templates, n_voxels_image);

          get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);
          // Rprintf("Got neighs.\n");

          // Rprintf("N[%u] = %u\n", 0, temp_neighs[0]);
          // Rprintf("N[%u] = %u\n", n_neighbours_patch - 1, temp_neighs[n_neighbours_patch - 1]);

          get_image_value(template4D.begin(), temp_neighs, temp_values, n_neighbours_patch);
          // Rprintf("Got image values.\n");

          double matchSum1 = 0;
          double matchSSQ1 = 0;

          double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);
          // double match = SSDPatch(input_values, temp_values, n_neighbours_patch, previous);
          //
          // Rprintf("Patch similarity.\n");

          if (match < similarities[k1 * actual_voxels + idx_voxel]) {

            // Rprintf("Improvement.\n");

            similarities[k1 * actual_voxels + idx_voxel] = match;
            kANN[k1 * actual_voxels + idx_voxel] = candidate;

          }

        }

        free(input_values);
        free(temp_values);
        free(input_patch);
        free(temp_neighs);

      }

    }

    // Rprintf("Dividing by 2 size.\n");

    search_size = search_size / 2;

    // free(search_neighbours);

  }

  // Rprintf("Exiting.\n");

  // free(result);

}


// [[Rcpp::export]]
void label_fusion_omp(IntegerVector labels4D,
                      int actual_voxels,
                      IntegerVector voxel_lookup_table,
                      IntegerVector label_ids,
                      IntegerVector kANN,
                      IntegerVector patch_neighbours,
                      int k,
                      double lambda,
                      double sigma2,
                      NumericVector match,
                      NumericVector new_voting,
                      int ncores = 2) {

  int n_neighbours_patch = patch_neighbours.size();
  int n_voxels_image = voxel_lookup_table.size();
  IntegerVector dims = labels4D.attr("dim");
  int n_labels = label_ids.size();

  double* similarities;
  int* labels;
  double* distances;

  int* patch_seg;

  int* input_patch;

  int* temp_neighs;

  int* counts = (int*) malloc(n_voxels_image * sizeof(int));
  memset(counts, 0, n_voxels_image * sizeof(int));

  int idx;

  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    idx = voxel_lookup_table[voxel];

    if (idx >= 0) {

      similarities = (double*) malloc(n_neighbours_patch * k * sizeof(double));
      labels = (int*) malloc(n_neighbours_patch * k * sizeof(int));
      distances = (double*) malloc(n_neighbours_patch * k * sizeof(double));

      input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));

      // Compute the graylevel of the input image in the patch around given voxel
      get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);

      // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(temp_neighs, patch_seg) schedule(static)
      for (int k1 = 0; k1 < k; k1++) {

        int best_match = kANN[k1 * actual_voxels + idx];
        int relative_index = best_match % n_voxels_image;

        double distance = 0.0;

        if (sigma2 > 0) {

          int x_voxel = voxel % dims[0];
          int x_candidate = relative_index % dims[0];
          int y_voxel = ((voxel - x_voxel) / dims[0]) % dims[1];
          int y_candidate = ((relative_index - x_candidate) / dims[0]) % dims[1];
          int z_voxel = ((voxel - x_voxel) / dims[0] - y_voxel) / dims[1];
          int z_candidate = ((relative_index - x_candidate) / dims[0] - y_candidate) / dims[1];

          distance = (x_voxel - x_candidate) * (x_voxel - x_candidate);
          distance += (y_voxel - y_candidate) * (y_voxel - y_candidate);
          distance += (z_voxel - z_candidate) * (z_voxel - z_candidate);

        }

        temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
        patch_seg = (int*) malloc(n_neighbours_patch * sizeof(int));

        get_neighbours(kANN[k1 * actual_voxels + idx], patch_neighbours.begin(), temp_neighs, n_neighbours_patch);
        get_label_value(labels4D.begin(), temp_neighs, patch_seg, n_neighbours_patch);

        double simil = match[k1 * actual_voxels + idx];

        // If we use the distances to neighboring patches in the search to weight the similarity
        for (int neigh = 0; neigh < n_neighbours_patch; neigh++) {

          similarities[k1 * n_neighbours_patch + neigh] = simil;
          labels[k1 * n_neighbours_patch + neigh] = patch_seg[neigh];
          distances[k1 * n_neighbours_patch + neigh] = distance;

        }

        free(temp_neighs);
        free(patch_seg);

      } // For all k


      // Normalize the similarities
      // First, compute the minimum similarity
      double min_simil = 1.e100;
      for (int temp = 0; temp < k; temp++) {

        if (similarities[temp * n_neighbours_patch] < min_simil) {

          min_simil = similarities[temp * n_neighbours_patch];

        }

      }
      // double min_simil = find_minimum(similarities, k * n_neighbours_patch);


      double h2 = lambda * min_simil + 1.e-15;
      for (int i = 0; i < k * n_neighbours_patch; i++) {

        similarities[i] = exp(-similarities[i] / h2 -  distances[i] / (sigma2 + 1.e-15));

      }

      // Compute the voting for each neighbour
      double cumulative_similarity = 0;
      for (int temp = 0; temp < k; temp++) {

        cumulative_similarity += similarities[temp * n_neighbours_patch];

      }

      // For all labels except the first (probably label == 0)
#pragma omp parallel for num_threads(ncores) schedule(dynamic)
      for (int lb = 1; lb < n_labels; lb++) {

        for (int i = 0; i < n_neighbours_patch; i++) {

          int which_voxel = input_patch[i];

          if (which_voxel > -1) {

            counts[which_voxel]++;

            double similarity = 0;

            for (int temp = 0; temp < k; temp++) {

              if (labels[temp * n_neighbours_patch + i] == label_ids[lb]) {

                similarity += similarities[temp * n_neighbours_patch + i];

              }

            }

            new_voting[lb * n_voxels_image + which_voxel] += similarity / cumulative_similarity;

          }

        }

      }

      free(input_patch);
      free(similarities);
      free(labels);
      free(distances);

    }

  }

#pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    double simil_0 = 1;
    for (int lb = 1; lb < n_labels; lb++) {

      simil_0 -= new_voting[lb * n_voxels_image + voxel];

    }

    new_voting[voxel] = simil_0;

  }

#pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    if (counts[voxel] > 0) {

      for (int lb = 0; lb < n_labels; lb++) {

        new_voting[lb * n_voxels_image + voxel] /= counts[voxel];

      }

    }

  }

  free(counts);

}


// [[Rcpp::export]]
void label_fusion2_omp(IntegerVector labels4D,
                       int actual_voxels,
                       IntegerVector voxel_lookup_table,
                       IntegerVector label_ids,
                       IntegerVector kANN,
                       IntegerVector patch_neighbours,
                       int k,
                       double lambda,
                       double sigma2,
                       NumericVector match,
                       NumericVector new_voting,
                       int ncores = 2) {

  int n_neighbours_patch = patch_neighbours.size();
  int n_voxels_image = voxel_lookup_table.size();
  IntegerVector dims = labels4D.attr("dim");
  int n_labels = label_ids.size();

  // double* similarities;
  // int* labels;
  // int* counts = (int*) malloc(n_voxels_image * sizeof(int));
  // memset(counts, 0, n_voxels_image * sizeof(int));
  int counts;
  int max_label = max(labels4D);

  int* new_labels = (int*) malloc((max_label + 1) * sizeof(int));

  for (int lb = 0; lb <= max_label; lb++) {

    new_labels[lb] = -1;

  }

  for (int lb = 0; lb < n_labels; lb++) {

    new_labels[label_ids[lb]] = lb;

  }

#pragma omp parallel for num_threads(ncores)
  for (int idx = 0; idx < actual_voxels; idx++) {

    double min_simil = 1.e10;
    for (int k1 = 0; k1 < k; k1++) {

      if (min_simil > match[k1 * actual_voxels + idx])
        min_simil = match[k1 * actual_voxels + idx];

    }

    double h2 = min_simil * lambda + 1.e-15;
    for (int k1 = 0; k1 < k; k1++) {

      match[k1 * actual_voxels + idx] = exp( match[k1 * actual_voxels + idx] / h2 );

    }

    double cumulative = 0;
    for (int k1 = 0; k1 < k; k1++) {

      cumulative += match[k1 * actual_voxels + idx];

    }

    for (int k1 = 0; k1 < k; k1++) {

      match[k1 * actual_voxels + idx] /= cumulative;

    }

  }


#pragma omp parallel for num_threads(ncores) private(counts) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    counts = 0;

    // Loop over all neighbours
    for (int neighbour = 0; neighbour < n_neighbours_patch; neighbour++) {

      int my_neighbour = voxel + patch_neighbours[neighbour];

      if ((my_neighbour < 0) || (my_neighbour >= n_voxels_image)) continue;

      // If idx_patch_center >= 0 is because it is the center of a patch
      int idx_patch_center = voxel_lookup_table[my_neighbour];

      // Rprintf("input_patch[%d] = %d, idx_patch_center = %u\n", neighbour, input_patch[neighbour], idx_patch_center);

      if (idx_patch_center >= 0) {

        // similarities = (double*) malloc(k * sizeof(double));
        // labels = (int*) malloc(k * sizeof(int));

        // memset(similarities, 0, k * sizeof(double));
        // memset(labels, 0, k * sizeof(int));

        int relative_difference = voxel - my_neighbour;

        for (int k1 = 0; k1 < k; k1++) {

          int candidate = kANN[k1 * actual_voxels + idx_patch_center];
          int voxel_candidate = candidate + relative_difference;

          // labels[k1] = labels4D[voxel_candidate];

          int lb = new_labels[labels4D[voxel_candidate]];

          if (lb >= 0)
            new_voting[lb * n_voxels_image + voxel] += match[k1 * actual_voxels + idx_patch_center];
          // for (int lb = 1; lb < n_labels; lb++) {
          //
          //   if (labels4D[voxel_candidate] == label_ids[lb]) {
          //
          //     // new_voting[lb * n_voxels_image + voxel] += similarities[k1] / cumulative_similarity;
          //
          //     new_voting[lb * n_voxels_image + voxel] += match[k1 * actual_voxels + idx_patch_center];
          //
          //     break;
          //
          //   }
          //
          // }

        }


        counts++;

        // free(similarities);
        // free(labels);


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


  // free(counts);

  free(new_labels);

}


// [[Rcpp::export]]
void label_fusion3_omp(IntegerVector labels4D,
                       int actual_voxels,
                       IntegerVector voxel_lookup_table,
                       IntegerVector label_ids,
                       IntegerVector kANN,
                       IntegerVector patch_neighbours,
                       int k,
                       double lambda,
                       double sigma2,
                       NumericVector match,
                       NumericVector new_voting,
                       int ncores = 2) {

  int n_neighbours_patch = patch_neighbours.size();
  int n_voxels_image = voxel_lookup_table.size();
  IntegerVector dims = labels4D.attr("dim");
  int n_labels = label_ids.size();

  double* similarities;
  int* labels;
  double* distances;

  int* patch_seg;

  int* input_patch;

  int* temp_neighs;
  int patch_size = (int)pow(n_neighbours_patch, 1/3);

  int* counts = (int*) malloc(n_voxels_image * sizeof(int));
  memset(counts, 0, n_voxels_image * sizeof(int));

  int idx;

  // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(temp_neighs, patch_seg, similarities, labels, distances, input_patch) schedule(dynamic)
  for (int init = 0; init < patch_size; init++) {

    for (int x = init; x < dims[0]; x+=patch_size) {

      for (int y = init; y < dims[1]; y+=patch_size) {

        for (int z = init; z < dims[2]; z+=patch_size) {

          int voxel = x + dims[0] * y + dims[0] * dims[1] * z;

          idx = voxel_lookup_table[voxel];

          if (idx >= 0) {

            similarities = (double*) malloc(n_neighbours_patch * k * sizeof(double));
            labels = (int*) malloc(n_neighbours_patch * k * sizeof(int));
            distances = (double*) malloc(n_neighbours_patch * k * sizeof(double));

            input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));

            // Compute the graylevel of the input image in the patch around given voxel
            get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);

            for (int k1 = 0; k1 < k; k1++) {

              int best_match = kANN[k1 * actual_voxels + idx];
              int relative_index = best_match % n_voxels_image;

              double distance = 0.0;

              if (sigma2 > 0) {

                int x_voxel = voxel % dims[0];
                int x_candidate = relative_index % dims[0];
                int y_voxel = ((voxel - x_voxel) / dims[0]) % dims[1];
                int y_candidate = ((relative_index - x_candidate) / dims[0]) % dims[1];
                int z_voxel = ((voxel - x_voxel) / dims[0] - y_voxel) / dims[1];
                int z_candidate = ((relative_index - x_candidate) / dims[0] - y_candidate) / dims[1];

                distance = (x_voxel - x_candidate) * (x_voxel - x_candidate);
                distance += (y_voxel - y_candidate) * (y_voxel - y_candidate);
                distance += (z_voxel - z_candidate) * (z_voxel - z_candidate);

              }

              temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
              patch_seg = (int*) malloc(n_neighbours_patch * sizeof(int));

              get_neighbours(kANN[k1 * actual_voxels + idx], patch_neighbours.begin(), temp_neighs, n_neighbours_patch);
              get_label_value(labels4D.begin(), temp_neighs, patch_seg, n_neighbours_patch);

              double simil = match[k1 * actual_voxels + idx];

              // If we use the distances to neighboring patches in the search to weight the similarity
              for (int neigh = 0; neigh < n_neighbours_patch; neigh++) {

                similarities[k1 * n_neighbours_patch + neigh] = simil;
                labels[k1 * n_neighbours_patch + neigh] = patch_seg[neigh];
                distances[k1 * n_neighbours_patch + neigh] = distance;

              }

              free(temp_neighs);
              free(patch_seg);

            } // For all k


            // Normalize the similarities
            // First, compute the minimum similarity
            double min_simil = 1.e100;
            for (int temp = 0; temp < k; temp++) {

              if (similarities[temp * n_neighbours_patch] < min_simil) {

                min_simil = similarities[temp * n_neighbours_patch];

              }

            }
            // double min_simil = find_minimum(similarities, k * n_neighbours_patch);


            double h2 = lambda * min_simil + 1.e-15;
            for (int i = 0; i < k * n_neighbours_patch; i++) {

              similarities[i] = exp(-similarities[i] / h2 -  distances[i] / (sigma2 + 1.e-15));

            }

            // Compute the voting for each neighbour
            double cumulative_similarity = 0;
            for (int temp = 0; temp < k; temp++) {

              cumulative_similarity += similarities[temp * n_neighbours_patch];

            }

            // For all labels except the first (probably label == 0)
            // #pragma omp parallel for num_threads(ncores) schedule(dynamic)
            for (int lb = 1; lb < n_labels; lb++) {

              for (int i = 0; i < n_neighbours_patch; i++) {

                int which_voxel = input_patch[i];

                if (which_voxel > -1) {

                  counts[which_voxel]++;

                  double similarity = 0;

                  for (int temp = 0; temp < k; temp++) {

                    if (labels[temp * n_neighbours_patch + i] == label_ids[lb]) {

                      similarity += similarities[temp * n_neighbours_patch + i];

                    }

                  }

                  new_voting[lb * n_voxels_image + which_voxel] += similarity / cumulative_similarity;

                }

              }

            }

            free(input_patch);
            free(similarities);
            free(labels);
            free(distances);

          }


        }

      }

    }

  }

#pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    double simil_0 = 1;
    for (int lb = 1; lb < n_labels; lb++) {

      simil_0 -= new_voting[lb * n_voxels_image + voxel];

    }

    new_voting[voxel] = simil_0;

  }

#pragma omp parallel for num_threads(ncores) schedule(dynamic, 10000)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    if (counts[voxel] > 0) {

      for (int lb = 0; lb < n_labels; lb++) {

        new_voting[lb * n_voxels_image + voxel] /= counts[voxel];

      }

    }

  }

  free(counts);

}

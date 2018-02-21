#include <Rcpp.h>
#include "nifti.h"
#include "similarity_measures.h"
#include "utils.h"
#include "neighbours.h"
#include <Rinternals.h>

// Add a flag to enable OpenMP at compile time
// 1[[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
// #ifdef _OPENMP
// #include <omp.h>
// #endif

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

using namespace Rcpp;

//##%######################################################%##
//#                                                          #
//####                     NIFTI I/O                      ####
//#                                                          #
//##%######################################################%##

// [[Rcpp::export]]
NumericVector fast_read_nifti(SEXP filename) {

  RNifti::NiftiImage image(filename);

  IntegerVector dims(image->ndim);
  NumericVector result(image->nvox);

  Rprintf("%u\n", image->nvox);
  Rprintf("%u\n", image->nbyper);
  Rprintf("%u\n", image->datatype);

  for (int i = 0; i < (int)image->nvox; i++) {

    result[i] = ((double*)image->data)[i];

  }
  // memcpy((void*)result.begin(), image->data, image->nbyper * image->nvox);
  // result.begin() = image->data;

  for (int i = 0; i < image->ndim; i++) {

    dims[i] = image->dim[i + 1];

  }

  result.attr("dim") = dims;

  return result;

}

template <typename T>
std::vector<T> nifti2array(RNifti::NiftiImage image) {

  return image.getData<T>();

}

template <int RTYPE, typename T>
Vector<RTYPE> nifti_list_to_array(std::vector< std::string > filenames,
                                  IntegerVector dims,
                                  std::vector<T> initial) {

  int nv = 1;
  for (int d = 0; d < dims.size(); d++) {

    nv *= dims[d];

  }
  int n1 = nv / filenames.size();

  Vector<RTYPE> res(nv);
  res.attr("dim") = dims;

  RNifti::NiftiImage image;

  for (int idx = 0; idx < n1; idx++) {

    res[idx] = initial[idx];

  }


  int last_index = n1;

  for (int i = 1; i < (int)filenames.size(); i++) {

    RNifti::NiftiImage image(filenames[i]);

    std::vector<T> tmp = image.getData<T>();
    for (int idx = last_index; idx < last_index + n1; idx++) {

      res[idx] = tmp[idx - last_index];

    }

    last_index += n1;

  }

  return res;

}


// [[Rcpp::export]]
SEXP read_nifti_batch_4d(std::vector< std::string > filenames) {

  RNifti::NiftiImage image(filenames[0]);

  IntegerVector dims(image->ndim + 1);

  SEXP result;

  for (int i = 0; i < image->ndim; i++) {

    dims[i] = image->dim[i + 1];

  }

  dims[image->ndim] = filenames.size();


  if ((image->datatype < 16) & (image->datatype > 1)) {

    result = nifti_list_to_array<INTSXP, int>(filenames, dims, image.getData<int>());

  } else {

    result = nifti_list_to_array<REALSXP, double>(filenames, dims, image.getData<double>());

  }

  return result;

}

// [[Rcpp::export]]
SEXP read_nifti(SEXP filename) {

  RNifti::NiftiImage image(filename);

  IntegerVector dims(image->ndim);

  SEXP result;

  for (int i = 0; i < image->ndim; i++) {

    dims[i] = image->dim[i + 1];

  }

  // result = image.toArray();

  if ((image->datatype < 16) & (image->datatype > 1)) {

    result = wrap(image.getData<int>());

  } else {

    result = wrap(image.getData<double>());

  }

  ((RObject)result).attr("dim") = dims;

  return result;

}



//##%######################################################%##
//#                                                          #
//####                        MALF                        ####
//#                                                          #
//##%######################################################%##

// [[Rcpp::export]]
void constrained_init_memory(IntegerVector dims,
                             int n_templates,
                             int patch_size,
                             int search_size,
                             IntegerVector voxel_lookup_table,
                             IntegerVector kANN,
                             int ncores = 1) {

  int actual_voxels;
  // IntegerVector dims = input_image.attr("dim");
  int n_voxels_image = dims[0] * dims[1] * dims[2]; //input_image.size();


  // Initialize similarities etc
  // First, indices for patch neighbours and search neighbours
  if (patch_size % 2 == 0) patch_size++;

  int k = n_templates;

  if (search_size % 2 == 0) search_size++;

  int limits = (int)((patch_size + search_size) / 2 + 1);

  Rprintf("n_voxels = %u, dims = [%u, %u, %u]\n", n_voxels_image, dims[0], dims[1], dims[2]);

  // Added an omp pragma directive to parallelize the loop with ncores
  // #pragma omp parallel for num_threads(ncores)
  for (int x = 0; x < dims[0]; x++) {

    for (int y = 0; y < dims[1]; y++) {

      for (int z = 0; z < dims[2]; z++) {

        int i = x + dims[0] * y + dims[0] * dims[1] * z;
        int idx = voxel_lookup_table[i];

        if (idx >= 0) {

          // Rprintf("x = %u, y = %u, z = %u, i = %u, idx = %d\n", x, y, z, i, idx);

          // Rprintf("i = %u, idx = %d, actual_voxels = %u\n", i, idx, actual_voxels);

          for (int k1 = 0; k1 < k; k1++) {

            // Rprintf("k1 = %u, k = %u\n", k1, k);

            int candidate = i;
            int template_id = k1 * 0; //rand_in_range(1, n_templates) - 1;
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

            // Rprintf("candidate = %u\n", candidate);
            // Rprintf("index_in_kANN = %u, length(kANN) = %u\n", k1 * actual_voxels + idx, n_templates * actual_voxels);

            kANN[k1 * actual_voxels + idx] = template_id * n_voxels_image + candidate;
            // kANN[k1 * actual_voxels + idx] = candidate;

            // Rprintf("Assigned\n");

          }

        }


      }

    }

  }

  // Rprintf("Exiting\n");

}


// [[Rcpp::export]]
void patches_similarity_memory(NumericVector input_image,
                               StringVector template_filenames,
                               int actual_voxels,
                               IntegerVector voxel_lookup_table,
                               IntegerVector patch_neighbours,
                               IntegerVector kANN,
                               NumericVector similarities,
                               int ncores = 1) {


  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  int k = template_filenames.size();
  // IntegerVector template_dims = template4D.attr("dim");

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  // int CHUNK_SIZE = n_voxels_image / ncores;

  for (int k1 = 0; k1 < k; k1++) {

    // Read actual image
    RNifti::NiftiImage image(template_filenames[k1]);

    std::vector<double> data = image.getData<double>();

    // Added an omp pragma directive to parallelize the loop with ncores
#pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(guided) shared(data)
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

        int candidate = kANN[k1 * actual_voxels + idx];
        Rprintf("k1 = %u, n_voxels_image = %u, voxel = %u, Candidate = %u\n", k1, n_voxels_image, voxel, candidate);

        get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

        for (int n = 0; n < n_neighbours_patch; n++) {

          temp_values[n] = data[temp_neighs[n]];

        }

        double matchSum1 = 0, matchSSQ1 = 0;

        double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);

        similarities[k1 * actual_voxels + idx] = match;

        free(input_values);
        free(temp_values);
        free(input_patch);
        free(temp_neighs);

      }

    }

  }

}


// [[Rcpp::export]]
void propagation_step_memory(NumericVector input_image,
                             StringVector template_filenames,
                             int actual_voxels,
                             IntegerVector voxel_lookup_table,
                             IntegerVector patch_neighbours,
                             IntegerVector kANN,
                             int direction,
                             int patch_size,
                             int stride,
                             NumericVector similarities,
                             int ncores = 1) {

  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");
  // IntegerVector template_dims = template4D.attr("dim");
  // int n_templates = template_dims[3];

  int k = template_filenames.size();

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int* neighbour_offset = (int*) malloc(3 * sizeof(int));
  neighbour_offset[0] = direction * stride;
  neighbour_offset[1] = neighbour_offset[0] * dims[0];
  neighbour_offset[2] = neighbour_offset[1] * dims[1];

  // int changes = 0;

  // Loop over all k ANN.
  for (int k1 = 0; k1 < k; k1++) {

    // Read actual image
    RNifti::NiftiImage image(template_filenames[k1]);

    std::vector<double> data = image.getData<double>();

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

        // Loop over all 6 neighbours
        for (int neighbour = 0; neighbour < 3; neighbour++) {

          int loc = voxel + neighbour_offset[neighbour];

          if (loc < 0 || loc >= n_voxels_image) continue;

          int idx_loc = voxel_lookup_table[loc];

          if (idx_loc >= 0) {

            int original_candidate = kANN[idx_loc + k1 * actual_voxels];

            int candidate = original_candidate - neighbour_offset[neighbour];

            int x_orig = candidate % dims[0];
            int y_orig = ((candidate - x_orig) / dims[0]) % dims[1];
            int z_orig = ((candidate - x_orig) / dims[0] - y_orig) / dims[1];

            if (x_orig < patch_size || x_orig > dims[0] - patch_size ||
                y_orig < patch_size || y_orig > dims[1] - patch_size ||
                z_orig < patch_size || z_orig > dims[2] - patch_size)
              continue;

            get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

            for (int n = 0; n < n_neighbours_patch; n++) {

              temp_values[n] = data[temp_neighs[n]];

            }

            double matchSum1 = 0, matchSSQ1 = 0;

            double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);

            if (match < similarities[k1 * actual_voxels + idx]) {

              similarities[k1 * actual_voxels + idx] = match;
              kANN[k1 * actual_voxels + idx] = candidate;

            }

          }

        }

        free(input_values);
        free(temp_values);
        free(input_patch);
        free(temp_neighs);

      }

    }

  }

  free(neighbour_offset);

}


// [[Rcpp::export]]
void constrained_random_search_memory(NumericVector input_image,
                                      StringVector template_filenames,
                                      int actual_voxels,
                                      IntegerVector voxel_lookup_table,
                                      IntegerVector kANN,
                                      int patch_size,
                                      IntegerVector patch_neighbours,
                                      int search_size_max,
                                      NumericVector similarities,
                                      int max_random_neighbours,
                                      int ncores = 1) {


  int n_voxels_image = input_image.size();
  IntegerVector dims = input_image.attr("dim");

  int k = template_filenames.size();

  int n_neighbours_patch = patch_neighbours.size();

  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  int search_size = search_size_max / 2;

  while (search_size > 1)  {

    // Rprintf("search_size = %u\n", search_size);

    if (search_size % 2 == 0) search_size++;

    // Loop over all k ANN.
    for (int k1 = 0; k1 < k; k1++) {

      // Read actual image
      RNifti::NiftiImage image(template_filenames[k1]);

      std::vector<double> data = image.getData<double>();

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

          // for (int k1 = 0; k1 < k; k1++) {

          int idx = kANN[idx_voxel + k1 * actual_voxels];

          int relative_index = idx;

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

          for (int n = 0; n < n_neighbours_patch; n++) {

            temp_values[n] = data[temp_neighs[n]];

          }
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

          // }

          free(input_values);
          free(temp_values);
          free(input_patch);
          free(temp_neighs);

        }

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
void label_fusion_memory(StringVector label_filenames,
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
  // IntegerVector dims = labels4D.attr("dim");
  int n_labels = label_ids.size();
  int k = label_filenames.size();

  // double* similarities;
  // int* labels;
  int* counts = (int*) malloc(n_voxels_image * sizeof(int));
  memset(counts, 0, n_voxels_image * sizeof(int));
  // int counts;
  int max_label = 10000; //max(labels4D);

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


  // Loop over all k ANN.
  for (int k1 = 0; k1 < k; k1++) {

    // Read actual image
    RNifti::NiftiImage image(label_filenames[k1]);

    std::vector<int> data = image.getData<int>();

#pragma omp parallel for num_threads(ncores) shared(counts) schedule(dynamic, 10000)
    for (int voxel = 0; voxel < n_voxels_image; voxel++) {

      // counts = 0;

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

          // for (int k1 = 0; k1 < k; k1++) {

          int candidate = kANN[k1 * actual_voxels + idx_patch_center];
          int voxel_candidate = candidate + relative_difference;

          // labels[k1] = labels4D[voxel_candidate];

          int lb = new_labels[data[voxel_candidate]];

          if (lb >= 0)
            new_voting[lb * n_voxels_image + voxel] += match[k1 * actual_voxels + idx_patch_center];

          // }


          counts[voxel]++;

          // free(similarities);
          // free(labels);


        } // Patch

      } // All neighbours




    } // All voxels


  }

#pragma omp parallel for num_threads(ncores)
  for (int voxel = 0; voxel < n_voxels_image; voxel++) {

    if (counts[voxel] > 0) {

      for (int lb = 1; lb < n_labels; lb++) {

        new_voting[lb * n_voxels_image + voxel] /= counts[voxel] * k;

      }

    }

    double simil_0 = 1;
    for (int lb = 1; lb < n_labels; lb++) {

      simil_0 -= new_voting[lb * n_voxels_image + voxel];

    }

    new_voting[voxel] = simil_0;

  }


  // free(counts);

  free(new_labels);

}

int count_eligible_voxels2(int ndims,
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


//[[Rcpp::export]]
List obtain_candidates_memory(NumericVector image,
                       StringVector template_files,
                       int patch_size,
                       int search_size,
                       int stride,
                       int max_iter,
                       int max_random_neighbours,
                       int ncores = 2) {

  IntegerVector dims = image.attr("dim");
  int n_voxels_image = image.size();

  IntegerVector voxel_lookup_table(n_voxels_image);

  int limits = (int)((patch_size + search_size) / 2 + 1);
  int ndims = dims.size();

  Rprintf("Init\n");

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

  int actual_voxels = count;

  Rprintf("Actual voxels = %u\n", actual_voxels);

  int k = template_files.size();

  IntegerVector kANN(k * actual_voxels);
  NumericVector similarities(k * actual_voxels);

  Rprintf("Size of kANN = %u\n", kANN.size());

  Rprintf("Constrained init\n");

  // Added an omp pragma directive to parallelize the loop with ncores
  for (int x = 0; x < dims[0]; x++) {

    for (int y = 0; y < dims[1]; y++) {

      for (int z = 0; z < dims[2]; z++) {

        int i = x + dims[0] * y + dims[0] * dims[1] * z;
        int idx = voxel_lookup_table[i];

        // if (idx >= 0) {
        //
        //   // Rprintf("x = %u, y = %u, z = %u, i = %u, idx = %d\n", x, y, z, i, idx);
        //
        //   // Rprintf("i = %u, idx = %d, actual_voxels = %u\n", i, idx, actual_voxels);
        //
        //   for (int k1 = 0; k1 < k; k1++) {
        //
        //     // Rprintf("k1 = %u, k = %u\n", k1, k);
        //
        //     int candidate = i;
        //     // int template_id = k1; //rand_in_range(1, n_templates) - 1;
        //     int x_displacement, y_displacement, z_displacement;
        //
        //     do {
        //
        //       x_displacement = rand_in_range(-search_size, search_size);
        //
        //     } while ((x + x_displacement <= limits) | (x + x_displacement > dims[0] - limits));
        //
        //     do {
        //
        //       y_displacement = rand_in_range(-search_size, search_size);
        //
        //     } while ((y + y_displacement <= limits) | (y + y_displacement > dims[1] - limits));
        //
        //     do {
        //
        //       z_displacement = rand_in_range(-search_size, search_size);
        //
        //     } while ((z + z_displacement <= limits) | (z + z_displacement > dims[2] - limits));
        //
        //     int absolute_displacement = x_displacement + dims[0] * y_displacement + dims[0] * dims[1] * z_displacement;
        //     candidate = i + absolute_displacement;
        //
        //     // Rprintf("candidate = %u\n", candidate);
        //     // Rprintf("index_in_kANN = %u, length(kANN) = %u\n", k1 * actual_voxels + idx, n_templates * actual_voxels);
        //
        //     // kANN[k1 * actual_voxels + idx] = template_id * n_voxels_image + candidate;
        //     kANN[k1 * actual_voxels + idx] = candidate;
        //
        //     // Rprintf("Assigned\n");
        //
        //   }
        //
        // }
        //

      }

    }

  }

  Rprintf("Neighbours\n");

  int n_neighbours_patch = pow(patch_size, 3);
  IntegerVector patch_neighbours(n_neighbours_patch);

  get_neighbours_indices(dims.begin(), ndims, patch_size, patch_neighbours.begin());

  Rprintf("Patches similarity\n");
  // Then, auxiliary variables
  double *input_values, *temp_values;
  int *input_patch, *temp_neighs;

  // int CHUNK_SIZE = n_voxels_image / ncores;

  for (int k1 = 0; k1 < k; k1++) {

    // Read actual image
    RNifti::NiftiImage nii_image(template_files[k1]);

    std::vector<double> data = nii_image.getData<double>();

    // Added an omp pragma directive to parallelize the loop with ncores
// #pragma omp parallel for num_threads(ncores) private(input_patch, input_values, temp_neighs, temp_values) schedule(guided) shared(data)
    for (int voxel = 0; voxel < n_voxels_image; voxel++) {

      int idx = voxel_lookup_table[voxel];

      if (idx >=0) {

        input_patch = (int*) malloc(n_neighbours_patch * sizeof(int));
        input_values = (double*) malloc(n_neighbours_patch * sizeof(double));

        temp_neighs = (int*) malloc(n_neighbours_patch * sizeof(int));
        temp_values = (double*) malloc(n_neighbours_patch * sizeof(double));

        get_neighbours(voxel, patch_neighbours.begin(), input_patch, n_neighbours_patch);
        get_image_value(image.begin(), input_patch, input_values, n_neighbours_patch);
        normalize(input_values, n_neighbours_patch);

        int candidate = kANN[k1 * actual_voxels + idx];
        // Rprintf("k1 = %u, n_voxels_image = %u, voxel = %u, Candidate = %u\n", k1, n_voxels_image, voxel, candidate);

        get_neighbours(candidate, patch_neighbours.begin(), temp_neighs, n_neighbours_patch);

        for (int n = 0; n < n_neighbours_patch; n++) {

          temp_values[n] = data[temp_neighs[n]];

        }

        double matchSum1 = 0, matchSSQ1 = 0;

        double match = patch_similarity(temp_values, input_values, n_neighbours_patch, matchSum1, matchSSQ1);

        similarities[k1 * actual_voxels + idx] = match;

        free(input_values);
        free(temp_values);
        free(input_patch);
        free(temp_neighs);

      }

    }

  }

  Rprintf("Start iterating.\n");
  // Rprintf("Mean similarity = ", mean(similarities), "\n");

  int direction = -1;

  for (int iter = 0; iter < max_iter; iter++) {

    propagation_step_memory(image,
                            template_files,
                            actual_voxels,
                            voxel_lookup_table,
                            patch_neighbours,
                            kANN,
                            direction,
                            patch_size,
                            stride,
                            similarities,
                            ncores);

    // cat("Mean similarity after PS = ", mean(similarities), "\n")

    direction *= -1;

    constrained_random_search_memory(image,
                                     template_files,
                                     actual_voxels,
                                     voxel_lookup_table,
                                     kANN,
                                     patch_size,
                                     patch_neighbours,
                                     search_size,
                                     similarities,
                                     max_random_neighbours,
                                     ncores);

    // cat("Mean similarity after CRS = ", mean(similarities), "\n")

  }

  return List::create(_["actual_voxels"] = actual_voxels,
                      _["voxel_lookup_table"] = voxel_lookup_table,
                      _["patch_neighbours"] = patch_neighbours,
                      _["candidates"] = kANN,
                      _["similarities"] = similarities);

}

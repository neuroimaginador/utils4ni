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


void which_max(double* input, int ndims, int* dims, double* output) {

  // We should have (ndims == kernel_ndims) or (ndims == kernel_ndims + 1)
  //

  if (ndims == 3) {

    for (int x = 0; x < dims[0]; x++) {

      for (int y = 0; y < dims[1]; y++) {

        double max;
        int pos = 0;

        for (int z = 0; z < dims[2]; z++) {

          int voxel = x + dims[0] * y + dims[0] * dims[1] * z;

          if (z == 0) {

            max = input[voxel];

          }

          if (max < input[voxel]) {

            max = input[voxel];
            pos = z;

          }

        }

        output[x + dims[0] * y] = (double)pos;

      }

    }

  } else {


    for (int x = 0; x < dims[0]; x++) {

      for (int y = 0; y < dims[1]; y++) {

        for (int z = 0; z < dims[2]; z++) {

          double max;
          int pos = 0;

          for (int k = 0; k < dims[3]; k++) {

            int voxel = x + dims[0] * y + dims[0] * dims[1] * z +  dims[0] * dims[1] * dims[2] * k;

            if (k == 0) {

              max = input[voxel];

            }

            if (max < input[voxel]) {

              max = input[voxel];
              pos = k;

            }

          }

          output[x + dims[0] * y + dims[0] * dims[1] * z ] = (double)pos;


        }

      }

    }

  }

}

//[[Rcpp::export]]
NumericVector which_max(NumericVector image) {

  IntegerVector dims = image.attr("dim");
  int ndims = dims.size();
  int image_size = image.size();
  int new_size = (int)(image_size / dims[ndims - 1]);
  IntegerVector new_dims(ndims - 1);
  for (int d = 0; d < ndims - 1; d++) {

    new_dims[d] = dims[d];

  }

  NumericVector segmentation(new_size);

  which_max(image.begin(), ndims, dims.begin(), segmentation.begin());

  segmentation.attr("dim") = new_dims;

  return segmentation;

}

void to_categorical_volume_cpp(double* image, int ndims, int* dims, int* unique_labels, int n_classes, int* segmentation) {

  if (ndims == 3) {

    for (int x = 0; x < dims[0]; x++) {

      for (int y = 0; y < dims[1]; y++) {

        for (int z = 0; z < dims[2]; z++) {

          int voxel = x + dims[0] * y + dims[0] * dims[1] * z;

          int k = (int)image[voxel];
          int target_voxel = x + dims[0] * y + dims[0] * dims[1] * z +  dims[0] * dims[1] * dims[2] * k;
          segmentation[target_voxel] = 1;

        }

      }

    }

  }

  if (ndims == 4) {

    // First dimension is batch_size

    for (int batch = 0; batch < dims[0]; batch++) {

      for (int x = 0; x < dims[1]; x++) {

        for (int y = 0; y < dims[2]; y++) {

          for (int z = 0; z < dims[3]; z++) {

            int voxel = batch + x * dims[0] + dims[0] * dims[1] * y + dims[0] * dims[1] * dims[2] * z;

            int k = (int)image[voxel];
            int target_voxel = voxel +  dims[0] * dims[1] * dims[2] * dims[3] * k;
            segmentation[target_voxel] = 1;

          }

        }

      }


    }

  }


}


//[[Rcpp::export]]
IntegerVector to_categorical_volume_cpp(NumericVector image, IntegerVector unique_labels) {

  IntegerVector dims = image.attr("dim");
  int ndims = dims.size();
  int image_size = image.size();

  int n_classes = unique_labels.size();

  IntegerVector new_dims(ndims + 1);
  int new_size = (int)(image_size * n_classes);

  for (int d = 0; d < ndims; d++) {

    new_dims[d] = dims[d];

  }
  new_dims[ndims] = n_classes;

  IntegerVector segmentation(new_size);

  to_categorical_volume_cpp(image.begin(), ndims, dims.begin(), unique_labels.begin(), n_classes, segmentation.begin());

  segmentation.attr("dim") = new_dims;

  return segmentation;


}

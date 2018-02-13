#include <Rcpp.h>

#include "RNifti.h"
#include "RNiftiAPI.h"
// #include "pointers.h"


using namespace Rcpp;

// SEXP rnifti_to_arraynd (SEXP image_) {
//
//   RNifti::NiftiImage image(image_);
//
//   struct array_nd *res = Calloc(1, struct array_nd);
//
//   res->n = image->nvox;
//   res->ndims = image->ndim;
//   res->datatype = image->nbyper == 4 ? DT_INT : DT_DOUBLE; // image->datatype;
//   res->data = image->data;
//   res->dims = (int*) malloc(res->ndims * sizeof(int));
//
//   for (int i = 0; i < res->ndims; i++) {
//     res->dims[i] = image->dim[i + 1];
//
//   }
//
//   SEXP ext_out = PROTECT(R_MakeExternalPtr(res, R_NilValue, R_NilValue));
//   R_RegisterCFinalizerEx(ext_out, _finalizer, TRUE);
//   UNPROTECT(1);
//
//   // Rprintf("%i, %f, %f, %f\n", foo->n, foo->data[0], foo->data[1], foo->data[2]);
//
//   return ext_out;
//
// }

// [[Rcpp::export]]
NumericVector fast_read_nifti(SEXP filename) {

  RNifti::NiftiImage image(filename);

  IntegerVector dims(image->ndim);
  NumericVector result(image->nvox);

  Rprintf("%u\n", image->nvox);
  Rprintf("%u\n", image->nbyper);

  for (int i = 0; i < image->nvox; i++) {

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

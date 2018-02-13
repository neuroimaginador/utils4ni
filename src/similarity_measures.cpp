#include <Rcpp.h>
using namespace Rcpp;
#include "similarity_measures.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double SSDPatch(double* PatchImg, double* PatchTemplate, int f) {

  /* SSD */
  double d = 0.0;

  for (int i = 0; i < f; i++) {

    d = d + (PatchImg[i] - PatchTemplate[i]) * (PatchImg[i] - PatchTemplate[i]);

  }

  d = d / f;

  return d;

}

double SSDPatch(double* PatchImg, double* PatchTemplate, int f, double previousSSD) {

  /* SSD */
  double d = 0.0;

  for (int i = 0; i < f; i++) {

    d = d + (PatchImg[i] - PatchTemplate[i]) * (PatchImg[i] - PatchTemplate[i]);

    if (d > previousSSD * f) {

      // Rprintf("Early exit %u/%u.\n", i, f);
      break;

    }

  }

  d = d / f;

  return d;

}

double patch_similarity(double* psearch, double* normtrg, int n, double &sum_psearch, double &ssq_psearch) {

  // Here the patch normtrg should already be normalized.
  // We simultaneously compute the patch stats and solve the problem
  double sum_uv = 0;
  for (int i = 0; i < n; i++) {
    double u = psearch[i];
    double v = normtrg[i];
    sum_psearch += u;
    ssq_psearch += u * u;
    sum_uv += u * v;
  }

  double var_u_unnorm = ssq_psearch - sum_psearch * sum_psearch / n;
  if(var_u_unnorm < 1.0e-6)
    var_u_unnorm = 1.0e-6;

  if(sum_uv > 0)
    return - (sum_uv * sum_uv) / var_u_unnorm;
  else
    return (sum_uv * sum_uv) / var_u_unnorm;

}


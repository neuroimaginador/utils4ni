#include <Rcpp.h>
#include <Rinternals.h>
using namespace Rcpp;

double SSDPatch(double* PatchImg, double* PatchTemplate, int f);
double SSDPatch(double* PatchImg, double* PatchTemplate, int f, double previousSSD);
double patch_similarity(double* psearch, double* normtrg, int n, double &sum_psearch, double &ssq_psearch);
double similarity(double* PatchImg, double* PatchTemplate, int n, int method);

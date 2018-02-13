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
NumericVector cvolfilter2d(NumericVector ima, IntegerVector dims, int v) {

  double media;
  int i,j,k,ii,jj,ni,nj,indice;

  // Result
  NumericVector fima(dims[0] * dims[1] * dims[2]);

  //precomputo de momentos de orden 1 y 2

  for (k = 0; k < dims[2]; k++) {

    for (i = 0; i < dims[1]; i++) {

      for (j = 0; j < dims[0]; j++) {

        if (ima[k * (dims[0] * dims[1]) + (i * dims[0]) + j]>0) {

          media = 0;
          indice = 0;

          for(ii = -v; ii <= v; ii++) {

            for (jj = -v; jj <= v; jj++) {

              ni = i + ii;
              nj = j + jj;

              if (ni >= 0 && nj >= 0 && ni < dims[1] && nj < dims[0]) {

                if (ima[k * (dims[0] * dims[1]) + (ni * dims[0]) + nj] > 0) {

                  media = media + ima[k * (dims[0] * dims[1]) + (ni * dims[0]) + nj];
                  indice = indice + 1;

                }

              }

            }

          }

          media = media / indice;
          fima[k * (dims[0] * dims[1]) + (i * dims[0]) + j] = media;

        }
        else fima[k * (dims[0] * dims[1]) + (i * dims[0]) + j] = 0;
      }
    }
  }

  fima.attr("dim") = dims;
  fima.attr("class") = "array";

  return fima;
}


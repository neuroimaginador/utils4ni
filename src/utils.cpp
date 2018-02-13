#include <Rcpp.h>
using namespace Rcpp;
#include "utils.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

void normalize(double* vector, int n) {

  if (n > 1) {

    double sum = 0; // Sum value

    // For loop, note cpp index shift to 0
    for (int i = 0; i < n; i++) {

      // Shorthand for sum = sum + x[i]
      sum += vector[i];

    }
    double mean = sum/n;

    sum = 0;
    for (int i = 0; i < n; i++) {

      sum += (vector[i] - mean) * (vector[i] - mean); // Square

    }
    double sd = sqrt(sum/(n-1));

    for (int i = 0; i < n; i++) {

      vector[i] = (vector[i] - mean) / (sd + 1.e-15);

    }
    // Rf_PrintValue(res);

  }

}

void print_timer(Timer timer) {

  Rf_PrintValue(timer);

  return;

}

int rand_in_range(int min, int max) {

  int value = min + (double)rand() / (double)RAND_MAX * (max - min + 1);

  return value;

}

//[[Rcpp::export]]
IntegerVector generate_random(int min, int max, int n) {

  IntegerVector res(n);

  for (int i = 0; i < n; i++) {

    res[i] = rand_in_range(min, max);

  }

  return res;

}

#include <Rcpp.h>
#include <Rdefines.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;

void normalize(double* vector, int n);
void print_timer(Timer timer);
int rand_in_range(int min, int max);

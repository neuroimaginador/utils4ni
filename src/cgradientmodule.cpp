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
NumericVector cgradientmodule(NumericVector ima, IntegerVector dims) {


  // Result
  NumericVector fima(dims[0] * dims[1] * dims[2]);
  double dx,dy,dz;
  int x1,x2,y1,y2,z1,z2,i,j,k;

  //precomputo de momentos de orden 1 y 2

  for(k=0;k<dims[2];k++)
  {
    for(i=0;i<dims[1];i++)
    {
      for(j=0;j<dims[0];j++)
      {
        if(ima[k*(dims[0]*dims[1])+(i*dims[0])+j]>0)
        {
          x1=((j+1)<dims[0])?(j+1):j;
          x2=(j-1)<0?0:(j-1);
          dx=(ima[k*(dims[0]*dims[1])+(i*dims[0])+x1]-ima[k*(dims[0]*dims[1])+(i*dims[0])+x2])/(x2-x1);

          y1=((i+1)<dims[1])?(i+1):i;
          y2=(i-1)<0?0:(i-1);
          dy=(ima[k*(dims[0]*dims[1])+(y1*dims[0])+j]-ima[k*(dims[0]*dims[1])+(y2*dims[0])+j])/(y2-y1);

          z1=((k+1)<dims[2])?(k+1):k;
          z2=(k-1)<0?0:(k-1);
          dz=(ima[z1*(dims[0]*dims[1])+(i*dims[0])+j]-ima[z2*(dims[0]*dims[1])+(i*dims[0])+j])/(z2-z1);

          fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=sqrt(dx*dx+dy*dy+dz*dz);
        }
        else fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=0;
      }
    }
  }

  fima.attr("dim") = dims;
  fima.attr("class") = "array";

  return fima;

}


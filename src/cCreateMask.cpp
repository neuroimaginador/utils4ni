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
NumericVector cCreateMask(NumericVector ima, IntegerVector dims, int c, NumericVector medias) {

  // Result
  NumericVector fima(dims[0] * dims[1] * dims[2]);

  //Declarations
  double mindis,dis;
  int i,j,k,ii,label;

  //precomputo de momentos de orden 1 y 2

  for(k=0;k<dims[2];k++)
  {
    for(i=0;i<dims[1];i++)
    {
      for(j=0;j<dims[0];j++)
      {
        if(ima[k*(dims[0]*dims[1])+(i*dims[0])+j]>0)
        {
          mindis=fabs(ima[k*(dims[0]*dims[1])+(i*dims[0])+j]-medias[0]);
          label=1;
          for(ii=1;ii<c;ii++)
          {
            dis=fabs(ima[k*(dims[0]*dims[1])+(i*dims[0])+j]-medias[ii]);
            if(dis<mindis)
            {
              mindis=dis;
              label=ii+1;
            }
          }

          fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=label;
        }
        else fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=0;
      }
    }
  }


  fima.attr("dim") = dims;
  fima.attr("class") = "array";

  return fima;

}


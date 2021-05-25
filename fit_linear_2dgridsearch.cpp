#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fit_linear_2dgridsearch(NumericMatrix x, NumericVector inter, float by) {
   NumericVector x1 = x(_, 0);
   NumericVector x2 = x(_, 1);
   
   unsigned int n = x1.length();
   unsigned int nt = floor(6.283185307179586231996 / by);
   
   float sum_inter = sum(inter);
   float theta, suml, sumr, best;
   
   NumericVector best_par(2);
   best_par[0] = 0.0;
   best_par[1] = std::numeric_limits<double>::max();
   best = sum_inter;
   
   NumericVector vals;
   std::vector<std::tuple<float,float>> zipped(n);

   for (int i = 0; i <= nt; i++) {
      theta = i * by;
      vals  = cos(theta) * x1 + sin(theta) * x2;
      
      for(int j = 0; j < n; j++) {
         std::get<0>(zipped[j]) = vals[j];
         std::get<1>(zipped[j]) = inter[j];
      }
      std::sort(zipped.begin(), zipped.end());
      
      suml = 0.0;
      sumr = sum_inter;

      for (int k = 0; k < n-1; k++) {
         suml += std::get<1>(zipped[k]);
         sumr -= std::get<1>(zipped[k]);
         if(sumr - suml > best){
            best = sumr - suml;
            best_par[0] = theta;
            best_par[1] = -(std::get<0>(zipped[k]) + std::get<0>(zipped[k+1])) / 2.0;
         }
      }
   }
   
   NumericVector out(3);
   out[0] = best_par[1];
   out[1] = cos(best_par[0]);
   out[2] = sin(best_par[0]);
   return out;
}
#include <Rcpp.h>
using namespace Rcpp;

// 计算 2^x - 1 的函数
// [[Rcpp::export]]
NumericVector exponential_base2_transform(NumericVector input) {
  NumericVector output(input.size());
  
  for(int i = 0; i < input.size(); i++) {
    output[i] = std::pow(2.0, input[i]) - 1.0;
  }
  
  return output;
}

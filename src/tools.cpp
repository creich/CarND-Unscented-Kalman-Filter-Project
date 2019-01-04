#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if( (estimations.size() == 0) || (ground_truth.size() == 0) ) {
     std::cerr << "ERROR: empty estimations OR ground_truth vector! .. quitting" << std::endl;
     return rmse;
  }

  //  * the estimation vector size should equal ground truth vector size
  if( estimations.size() != ground_truth.size() ) {
     std::cerr << "ERROR: estimations vector and ground_truth vector have different sizes! .. quitting" << std::endl;
     return rmse;
  }

  // accumulate squared residuals
  VectorXd squared_res(4);
  squared_res << 0, 0, 0, 0;

  for (int i=0; i < estimations.size(); ++i) {
     VectorXd c = estimations[i] - ground_truth[i];
     VectorXd s_res = c.array() * c.array();
     squared_res = squared_res + s_res;
  }

  // calculate the mean
  VectorXd mean = squared_res / estimations.size();

  // calculate the squared root
  rmse = mean.array().sqrt();

  // return the result
  return rmse;
}

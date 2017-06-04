#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	VectorXd rmse(4);
	rmse.fill(0);


	if (estimations.size() == 0
	   || estimations.size() != ground_truth.size() ){
	  return rmse;
	}

	//accumulate squared residuals

	for(unsigned int i=0; i < estimations.size(); ++i){
	  // ... your code here
	  VectorXd diff = (estimations[i] - ground_truth[i]);
	  VectorXd diff_sqrd = (diff.array() * diff.array());
	  rmse += diff_sqrd;
	}

	//calculate the mean
	// ... your code here
	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();

	//calculate the squared root
	// ... your code here

	//return the result
	return rmse;
}

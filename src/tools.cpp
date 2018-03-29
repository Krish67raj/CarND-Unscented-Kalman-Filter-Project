#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */VectorXd rmse(estimations[0].size());
	
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	
	//accumulate squared residuals
	
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd dif = estimations[i] - ground_truth[i];
		dif = dif.array()*dif.array();
		rmse += dif;
	}
	//calculate the mean
	rmse = rmse/estimations.size();
	cout<<"size estimation      "<<estimations.size()<<endl;
	cout<<rmse<<endl;
	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
  
  
}
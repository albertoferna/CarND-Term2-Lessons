#include <iostream>
#include "Dense"
#include <vector>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd CalculateJacobian(const VectorXd& x_state);

int main() {

	/*
	 * Compute the Jacobian Matrix
	 */

	//predicted state  example
	//px = 1, py = 2, vx = 0.2, vy = 0.4
	VectorXd x_predicted(4);
	x_predicted << 1, 2, 0.2, 0.4;

	MatrixXd Hj = CalculateJacobian(x_predicted);

	cout << "Hj:" << endl << Hj << endl;

	return 0;
}

MatrixXd CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE
	float p_numerator = px * px + py * py;

	//check division by zero
	if (fabs(p_numerator) < 0.00000001) {
	    cout << "CalculateJacobian() - Error - Division by Zero" << endl;
	    return Hj;
	}
	//compute the Jacobian matrix
	Hj(0,0) = Hj(2,2) = px / sqrt(p_numerator);
	Hj(0,1) = Hj(2,3) = py / sqrt(p_numerator);
	Hj(1,0) = -py / p_numerator;
	Hj(1,1) = px / p_numerator;
	Hj(2,0) = py * (vx * py - vy * px) / pow(p_numerator, 3/2);
	Hj(2,0) = px * (vy * px - vx * py) / pow(p_numerator, 3/2);
	return Hj;
}

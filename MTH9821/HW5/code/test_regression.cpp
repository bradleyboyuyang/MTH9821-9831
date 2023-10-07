#include<iostream>
#include"regression.h"
#include <Eigen/Core>

using arr = Eigen::VectorXd;


int main(){

	arr coefficients(4);
    coefficients << 1.0, -2.0, 3.0, -4.0; // Set your polynomial coefficients

    arr x_values(3);
    x_values << 1.0, 2.0, 3.0; // Set the points at which to evaluate the polynomial
    
    arr result = evaluate_polynomial(coefficients, x_values);
    std::cout << "Result: " << std::endl << result << std::endl;
    return 0;
}


#include<iostream>
#include"regression.h"
#include <Eigen/Core>

using arr = Eigen::VectorXd;
using mat = Eigen::MatrixXd;


void evaluate_polynomial(){
    arr coefficients(4);
    coefficients << 1.0, -2.0, 3.0, -4.0; // Set your polynomial coefficients

    arr x_values(3);
    x_values << 1.0, 2.0, 3.0; // Set the points at which to evaluate the polynomial
    
    arr result = linalg::regression::evaluate_polynomial(coefficients, x_values);
    std::cout << "Result: " << std::endl << result << std::endl;
}

void hermite_polynomial(){
    // Define test coefficients and x values
    arr coefficients(5);
    coefficients << 1, 2, 3, 4, 5; // Example coefficients for the Hermite polynomial
    arr x(5);
    x << -1, 0, 1, 2, 3; // Example x values

    // Call the function with test data
    arr result = linalg::regression::evaluate_hermite_polynomial(coefficients, x);

    // Output the result
    std::cout << "Input coefficients: " << coefficients.transpose() << std::endl;
    std::cout << "Input x values: " << x.transpose() << std::endl;
    std::cout << "Result: " << result.transpose() << std::endl;
}

void test_vandermonde() {
    arr x(3);
    x << 1, 2, 3;
    mat expected(3, 3);
    expected << 1, 1, 1,
                1, 2, 4,
                1, 3, 9;
    
    mat result = linalg::regression::vandermonde(x, 2);
    
    // Comparing if result is approximately equal to expected (you can adjust the precision as needed)
    if((result - expected).norm() < 1e-9){
        std::cout << "Vandermot test worked" << std::endl;
    }else{
        std::cout << "Vandermot test did not work" << std::endl;
    }

}

void test_evaluate_hermite_polynomial_standardized() {
    // Define test inputs
    arr coefficients(3);
    coefficients << 1, 2, 3;
    arr x(3);
    x << 1, 2, 3;
    double mu = 2;
    double sigma = 1;

    // Define expected outputs (you should calculate these based on your inputs)
    arr expected(3);
    // The expected values are calculated manually or using another reliable method
    // for demonstration purposes, let's assume the expected results for standardized
    // inputs -1, 0, 1 are -2, 1, and 6 respectively.
    expected << -2, 1, 6;

    // Call the function and get the result
    arr result = linalg::regression::evaluate_hermite_polynomial_standardized(coefficients, x, mu, sigma);

    // Compare result to expected output with some tolerance
    if((result - expected).norm() < 1e-6){
        std::cout << "hermit polynomial test worked" << std::endl;
    }else{ 
        std::cout << (result - expected).norm() << std::endl;
        std::cout << "hermit polynomial did not work" << std::endl;

    }
}

using arr2 = Eigen::MatrixXd;

bool test_fit_linear_regression() {
    // Define test inputs
    arr2 X(3, 2);
    X << 1, 1,
          2, 1,
          3, 1;

    arr y(3);
    y << 1, 2, 3;

    // Call the function with preconditioning set to false
    auto result = linalg::regression::fit_linear_regression(X, y, false);
    
    // Extract beta and fitted values from the result
    arr beta = result.first;
    arr fitted = result.second;

    // Define expected outputs
    arr expected_beta(2);
    expected_beta << 1, 0;  // Expecting a slope of 1 and intercept of 0 for this simple case

    arr expected_fitted(3);
    expected_fitted << 1, 2, 3;  // Expecting fitted values to match original y for this simple case

    // Tolerance for floating-point comparison
    double tol = 1e-9;

    // Compare result to expected output with some tolerance
    bool beta_close = (beta - expected_beta).cwiseAbs().maxCoeff() < tol;
    bool fitted_close = (fitted - expected_fitted).cwiseAbs().maxCoeff() < tol;

    if(beta_close && fitted_close){
        std::cout << "passed test" << std::endl;
    }else{
        std::cout << "did not pass test" << std::endl;
    }
}






int main(){

    evaluate_polynomial();
    hermite_polynomial();
    test_vandermonde();
    test_evaluate_hermite_polynomial_standardized();
    test_fit_linear_regression();
    return 0;
}


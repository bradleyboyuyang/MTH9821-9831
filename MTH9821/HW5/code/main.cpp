#include<iostream>
#include"regression.h"
#include <Eigen/Core>




void check_polynomial(){
    // Define a set of coefficients for a polynomial: 2x^2 + 3x + 1
    Eigen::VectorXd coefficients(3);
    coefficients << 1, 3, 2;

    // Define a set of x values
    Eigen::VectorXd x_values(3);
    x_values << -1, 0, 1;

    // Evaluate the polynomial for the given x values
    Eigen::VectorXd result = linalg::regression::evaluate_polynomial(coefficients, x_values);

    // Print the result
    std::cout << "Result: " << result.transpose() << std::endl;

    // Check the result
    // For x = -1, the polynomial value is 2*(-1)^2 + 3*(-1) + 1 = 0
    // For x = 0, the polynomial value is 2*(0)^2 + 3*(0) + 1 = 1
    // For x = 1, the polynomial value is 2*(1)^2 + 3*(1) + 1 = 6
    Eigen::VectorXd expected_result(3);
    expected_result << 0, 1, 6;
    if ((result - expected_result).norm() < 1e-6) {
        std::cout << "Polynomial Evaluation works!" << std::endl;
    } else {
        std::cout << "Test failed!" << std::endl;
        std::cout << "Expected: " << expected_result.transpose() << std::endl;
        std::cout << "Obtained: " << result.transpose() << std::endl;
    }
}

void hermite_polynomial(){
    // Define test coefficients and x values
    Eigen::VectorXd coefficients(5);
    coefficients << 1, 2, 3, 4, 5; // Example coefficients for the Hermite polynomial
    Eigen::VectorXd x(5);
    x << -1, 0, 1, 2, 3; // Example x values

    // Call the function with test data
    Eigen::VectorXd result = linalg::regression::evaluate_hermite_polynomial(coefficients, x);

    // Output the result
    std::cout << "Input coefficients: " << coefficients.transpose() << std::endl;
    std::cout << "Input x values: " << x.transpose() << std::endl;
    std::cout << "Result: " << result.transpose() << std::endl;
}

void test_vandermonde() {
    Eigen::VectorXd x(3);
    x << 1, 2, 3;
    Eigen::MatrixXd expected(3, 3);
    expected << 1, 1, 1,
                1, 2, 4,
                1, 3, 9;
    
    Eigen::MatrixXd result = linalg::regression::vandermonde(x, 2);
    
    // Comparing if result is approximately equal to expected (you can adjust the precision as needed)
    if((result - expected).norm() < 1e-9){
        std::cout << "Vandermot test worked" << std::endl;
    }else{
        std::cout << "Vandermot test did not work" << std::endl;
    }

}

void test_evaluate_hermite_polynomial_standardized() {
    // Define test inputs
    Eigen::VectorXd coefficients(3);
    coefficients << 1, 2, 3;
    Eigen::VectorXd x(3);
    x << 1, 2, 3;
    double mu = 2;
    double sigma = 1;

    // Define expected outputs (you should calculate these based on your inputs)
    Eigen::VectorXd expected(3);
    // The expected values are calculated manually or using another reliable method
    // for demonstration purposes, let's assume the expected results for standardized
    // inputs -1, 0, 1 are -2, 1, and 6 respectively.
    expected << -2, 1, 6;

    // Call the function and get the result
    Eigen::VectorXd result = linalg::regression::evaluate_hermite_polynomial_standardized(coefficients, x, mu, sigma);

    // Compare result to expected output with some tolerance
    if((result - expected).norm() < 1e-6){
        std::cout << "hermit polynomial test worked" << std::endl;
    }else{ 
        std::cout << (result - expected).norm() << std::endl;
        std::cout << "hermit polynomial did not work" << std::endl;

    }
}

void test_fit_linear_regression(){
    // Create example data
    Eigen::MatrixXd X(5, 2);
    X << 1.0, 2.0,
         2.0, 3.0,
         3.0, 4.0,
         4.0, 5.0,
         5.0, 6.0;
         
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    // Fit the linear regression model
    auto result = linalg::regression::fit_linear_regression(X, y, true);
    Eigen::VectorXd beta = result.first;
    Eigen::VectorXd fitted_values = result.second;

    // Print the results
    std::cout << "Beta coefficients: \n" << beta << std::endl;
    std::cout << "Fitted values: \n" << fitted_values << std::endl;
    //return std::make_pair(beta, fitted_values);
}





int main(){
    check_polynomial();
    hermite_polynomial();
    test_vandermonde();
    test_evaluate_hermite_polynomial_standardized();
    test_fit_linear_regression();
    return 0;
}


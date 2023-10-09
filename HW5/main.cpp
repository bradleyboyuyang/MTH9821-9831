#include "regression.h"
#include "mc_regression.h"
#include "american_pricers.h"

#include <iostream>
#include <random>

int main(){

	typedef Eigen::ArrayXd arr;
	typedef Eigen::ArrayXXd arr2;
    using namespace montecarlo;
    std::default_random_engine gen;
	arr inputarr(4);
	inputarr(0) = -1.0; inputarr(1) = 0.0; inputarr(2) = 1.0; inputarr(3) = 4.0;
	arr inputarr2(4);
	inputarr2(0) = 3.0; inputarr2(1) = 9.0; inputarr2(2) = -2.0; inputarr2(3) = 6.0;
	arr2 test = std::get<0>(linalg::regression::hermite_vandermonde_standardized(inputarr, 3));
    arr2 test2 = (linalg::regression::hermite_vandermonde(inputarr, 3));
    arr2 test3 = linalg::regression::vandermonde(inputarr, 3);
	auto newtest = linalg::regression::fit_linear_regression(test, inputarr2, true);
	std::cout << newtest.first << std::endl;
	std::cout << newtest.second << std::endl;
//	std::cout << test << std::endl;
//    std::cout << test2 << std::endl;
//    std::cout << test3 << std::endl;
    //std::cout << newtest << std::endl;
	//std::cout << linalg::regression::evaluate_polynomial(inputarr, inputarr2) << std::endl;
	//std::cout << linalg::regression::evaluate_polynomial2(inputarr, inputarr2);

    const double spot = 100;
    const double strike = 100;
    const double maturity = 1;
    const double interest_rate = 0.05;
    const int w = 1;

    const int N = 10000;
    const int M = 100;
    const double dt = maturity / M;
    const double sigma = 0.2;

    //generate paths
    arr2 paths(N, M );
    paths.col(0).fill(spot);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M-1; j++) {
            std::normal_distribution<double> z(0.0, 1.0);
            paths(i, j + 1) = paths(i, j) * exp((interest_rate - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * z(gen));
        }
    }
    montecarlo::regression::PolynomialMCRegression mc_regression(M, 3, true);

    double P = regression_pricer_backward(spot, paths, w, strike, maturity, interest_rate,
                                          mc_regression, MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    std::cout << "Tsitsiklis Van Roy: " << P << std::endl;
    double P2 = regression_pricer_backward(spot, paths, w, strike, maturity, interest_rate,
                                          mc_regression, MonteCarloRegressionMethod::Longstaff_Schwartz);
    std::cout << "Longstaff Schwartz: " << P2 << std::endl;

    double P3 = regression_pricer_forward(spot, paths, w, strike, maturity, interest_rate,
                                          mc_regression);
    std::cout << "Forward: " << P3 << std::endl;
}
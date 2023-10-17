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
	auto newtest = linalg::regression::fit_linear_regression(test3, inputarr2, true);
	std::cout << newtest.first << std::endl;
	std::cout << newtest.second << std::endl;
//	std::cout << test << std::endl;
//    std::cout << test2 << std::endl;
//    std::cout << test3 << std::endl;
    //std::cout << newtest << std::endl;
	//std::cout << linalg::regression::evaluate_polynomial(inputarr, inputarr2) << std::endl;
	//std::cout << linalg::regression::evaluate_polynomial2(inputarr, inputarr2);

    const double spot = 40.5;
    const double strike = 44;
    const double maturity = 0.6;
    const double interest_rate = 0.04;
    const int w = -1;


    const int K = 100;
    const int N = 10000;
    const int M = 6;
    const double dt = maturity / M;
    const double sigma = 0.2;

    arr P(K);
    arr P1(K);
    arr P2(K);
    arr P3(K);

    for(int k = 0; k < K; k++){
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

        P(k) = regression_pricer_backward(spot, paths, w, strike, maturity, interest_rate,
                                              mc_regression, MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
        P1(k) = regression_pricer_forward(spot, paths, w, strike, maturity, interest_rate,
                                              mc_regression);
        P2(k) = regression_pricer_backward(spot, paths, w, strike, maturity, interest_rate,
                                               mc_regression, MonteCarloRegressionMethod::Longstaff_Schwartz);
        P3(k) = regression_pricer_forward(spot, paths, w, strike, maturity, interest_rate,
                                              mc_regression);
    }
    double P_mean = P.mean();
    double P1_mean = P1.mean();
    double P2_mean = P2.mean();
    double P3_mean = P3.mean();

    std::cout << "Tsitsiklis Van Roy: " << P_mean << std::endl;
    std::cout << "Tsitsiklis Van Roy Forward: " << P1_mean << std::endl;
    std::cout << "Longstaff Schwartz: " << P2_mean << std::endl;
    std::cout << "Longstaff Schwartz Forward: " << P3_mean << std::endl;


}
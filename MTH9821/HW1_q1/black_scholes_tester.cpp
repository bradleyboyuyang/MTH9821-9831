#include<iostream>
#include"black_scholes_pricer.hpp"

int main() {
    const double S0 = 41.0;        // Spot price
    const double K = 42.0;         // Strike price
    const double T = 0.75;         // Time to maturity (in years)
    const double sigma = 0.25;     // Volatility
    const double q = 0.01;         // Dividend yield
    const double r = 0.03;         // Risk-free interest rate
    double callPrice = blackScholesCall(S0, K, T, sigma, q, r);
    double putPrice = blackScholesPut(S0, K, T, sigma, q, r);

    std::cout << "European Call Option Price: " << callPrice << std::endl;
    std::cout << "European Put Option Price: " << putPrice << std::endl;

    return 0;
}

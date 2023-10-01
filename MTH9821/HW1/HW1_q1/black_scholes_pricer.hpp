#ifndef BlackScholes
#define BlackScholes
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

// Function to calculate the cumulative distribution function of the standard normal distribution
double N(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Black-Scholes formula for European call option
double blackScholesCall(double S0, double K, double T, double sigma, double q, double r) {
    double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S0 * std::exp(-q * T) * N(d1) - K * std::exp(-r * T) * N(d2);
}

double DeltablackScholesCall(double S0, double K, double T, double sigma, double q, double r) {
    double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return std::exp(-q * T) * N(d1);
}

double VegablackScholesCall(double S0, double K, double T, double sigma, double q, double r) {
    double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return S0*std::exp(-q*T)*std::sqrt(T)/std::sqrt(2*M_PI) * std::exp(-d1*d1/2);
}

// Black-Scholes formula for European call option
double blackScholesPut(double S0, double K, double T, double sigma, double q, double r) {
    double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return -S0 * std::exp(-q * T) * N(-d1) + K * std::exp(-r * T) * N(-d2);
}

double DeltablackScholesPut(double S0, double K, double T, double sigma, double q, double r) {
    double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return -std::exp(-q * T) * N(-d1);
}

double VegablackScholesPut(double S0, double K, double T, double sigma, double q, double r) {
    double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return S0*std::exp(-q*T)*std::sqrt(T)/std::sqrt(2*M_PI) * std::exp(-d1*d1/2);
}


#endif

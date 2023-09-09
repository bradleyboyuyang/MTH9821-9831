#include<cmath>
#include<iostream>
#include <fstream>
#include<tuple>
#include"black_scholes_pricer.hpp"
#include"bsm.hpp"
#include"CallPricer.hpp"

int main(){

    const double S0 = 41.0;        // Spot price
    const double K = 42.0;         // Strike price
    const double T = 0.75;         // Time to maturity (in years)
    const double sigma = 0.25;     // Volatility
    const double q = 0.01;         // Dividend yield
    const double r = 0.03;         // Risk-free interest rate
    double callPrice = blackScholesCall(S0, K, T, sigma, q, r);
    double DeltaCall = DeltablackScholesCall(S0, K,  T, sigma, q, r);
    double VegaCall = VegablackScholesCall( S0, K,  T,  sigma,  q, r);

    //store to csv file
    std::ofstream outputFile("resultsHW1.csv");
    unsigned int N = 10000;
    for(int i =0; i < 10; i++){
        auto [S,C,Delta,Vega] = Call_pricer(N,S0,K,T,sigma,q,r);
        outputFile << N << ";" << C - callPrice << ";" <<Delta - DeltaCall << ";" << Vega - VegaCall <<  std::endl; 
        std::cout << N << " " << C - callPrice << " " <<Delta - DeltaCall << " " << Vega - VegaCall << " " <<  std::endl; 
        N*=2;
    }
    return 0;
}
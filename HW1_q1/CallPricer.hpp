#ifndef CallPricer
#define CallPricer
#include<cmath>
#include<iostream>
#include<tuple>
#include"black_scholes_pricer.hpp"
#include"bsm.hpp"

std::tuple<double, double, double,double> Call_pricer(unsigned int N, double S0, double K, double T, double sigma,  double q, double r){
    double SSum =0;
    double CSum =0;
    double DeltaSum = 0;
    double VegaSum = 0;
    for(int i = 1; i < N +2; i++){
        double zi = bsm(double(i)/(N + 2));
        double Si = S0* std::exp((r - q - sigma*sigma/2) * T + sigma * std::sqrt(T) *zi);
        
        double Ci = 0;
        double Di =0;
        double Vegai = 0;
        if(Si > K){
            Vegai = Si*std::exp(-r*T)*(- sigma * T + std::sqrt(T) * zi);
            Ci = std::exp(- r * T) * (Si - K);
            Di = std::exp(-r*T)*Si/S0;
        }
        SSum += Si;
        CSum += Ci;
        DeltaSum += Di;
        VegaSum +=Vegai;
    }

    SSum /= (N);
    CSum /= (N);
    DeltaSum /= (N);
    VegaSum /= (N);
    return {SSum, CSum, DeltaSum,VegaSum};
}
#endif

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
    for(int i = 1; i < N; i++){
        double zi = bsm(double(i)/N);
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

    SSum /= (N-2);
    CSum /= (N-2);
    DeltaSum /= (N-2);
    VegaSum /= (N-2);
    return {SSum, CSum, DeltaSum,VegaSum};
}
#endif
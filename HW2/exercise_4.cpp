#include<iostream>
#include<vector>
#include<cmath>
#include <random>

#include"black_scholes_pricer.hpp"

//Need to check the timesteps
//
double lambda =4;
double sqV = 0.35;
double eta = 0.25;
double rho = -0.15;
double T = 0.5;
double K = 50;
double r = 0.05;
double dt = T/(2* 87500*32);

double max(double a, double b){
	if ( a < b){
		return b;
	}
	return a;
}

int main(){

    // Create a random number generator engine
    std::random_device rd;
    std::mt19937 gen(rd());  // Mersenne Twister engine for randomness

    // Define the mean and standard deviation of the normal distribution
    double mean = 0.0;          // Mean
    double stddev = 1.0;        // Standard Deviation

    // Create a normal distribution
    std::normal_distribution<double> distribution(mean, stddev);
 
	unsigned int n = 1000;
	unsigned int N = 2* 87500*32;
	double Vbar = 0.35*0.35;

    std::vector<double> V(n);
    std::vector<double> S(n);
    std::vector<double> OpVal(n);
	//Given parameters	
	S[0] = 50;
   	V[0] = 0.09;
	OpVal[0]=0;	
	for(int i = 1; i < n; i++){
        double z1 = distribution(gen);
        double z2 = distribution(gen);
        //double Vbar = std::accumulate(V.begin(), V.begin() +i , 0.)/V.size();
		S[i] = S[i-1] * std::exp( (r - max(V[i-1], 0)/2)*dt + std::sqrt(max(V[i-1], 0) )* std::sqrt(dt)*z1);
        V[i] = max(V[i-1], 0) - lambda *( max(V[i-1], 0) -Vbar) * dt + eta * max(V[i-1], 0) * std::sqrt(dt) * ( rho* z1  + std::sqrt(1  - rho*rho) * z2);
		OpVal[i] = blackScholesCall(S[i], 50, 1./2., std::sqrt(V[i]), 0.0, 0.05);
        //std::cout << S[i] << " " << V[i] << " " << Vbar << std::endl;
    }
       
	std::cout << std::accumulate(OpVal.begin(), OpVal.end(), 0.0)/n << std::endl;	
    //double Vbar = std::accumulate(V.begin(), V.end(), 0)/(n - 1);
std::cout << Vbar << std::endl;
    return 0;
}
		

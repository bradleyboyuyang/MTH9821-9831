#include<iostream>
#include<vector>
#include<cmath>
#include<random>
#include"black_scholes_pricer.hpp"

double max(double a, double b){
	if ( a < b){
		return b;
	}
	return a;
}

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

// Linear Congruential Generator parameters
const unsigned long long a = 39373;
const unsigned long long c = 0;
const unsigned long long k = 2147483647; // 2^31 - 1

double generateUniform() {
    static unsigned long long xi = 1;
    xi = (a * xi + c) % k;
    return static_cast<double>(xi) / k;
}

void boxMuller(double& z1, double& z2) {
    double u1 = generateUniform();
    double u2 = generateUniform();

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;

    z1 = r * cos(theta);
    z2 = r * sin(theta);
}


int main(){
    double lambda = 4;
    double sqV = 0.35;
    double eta = 0.25;
    double rho = -0.15;
    double T = 0.5;
    double K = 50;
    double r = 0.05;

    // Create a random number generator engine
    std::random_device rd;
    std::mt19937 gen(rd());  // Mersenne Twister engine for randomness

    // Define the mean and standard deviation of the normal distribution
    double mean = 0.0;          // Mean
    double stddev = 1.0;        // Standard Deviation

    // Create a normal distribution
    std::normal_distribution<double> distribution(mean, stddev);
 
	double Vbar = 0.35*0.35;

    //TODO:change time frame
    for(int k = 0; k <=5; k++){
        int n = 500* std::pow(2, k);
        std::cout << n << std::endl;
        std::vector<double> V(n);
        std::vector<double> S(n);
        std::vector<double> stor;
        std::vector<double> vola;
        std::vector<double> under;
	    S[0] = 50;
   	    V[0] = 0.09;
        for(int j=0; j < n; j++){
            double T = 0.5;
            double dt = 0.5/175.;
            for(int i = 1; i < 87; i++){

                double z1, z2;
                boxMuller(z1, z2);
                //double z1 = distribution(gen);
                //double z2 = distribution(gen);
		        S[i] = S[i-1] * std::exp( (r - max(V[i-1], 0)/2)*dt + std::sqrt(max(V[i-1], 0) )* std::sqrt(dt)*z1);
                //std::cout << "S " << S[i] << std::endl;
                V[i] = max(V[i-1], 0) - lambda *( max(V[i-1], 0) -Vbar) * dt + eta * std::sqrt(max(V[i-1], 0)) * std::sqrt(dt) * ( rho* z1  + std::sqrt(1  - rho*rho) * z2);
                //std::cout << "V " << V[i] << std::endl;
                //std::cout << "blackscholes " << blackScholesCall(S[i], K,  T, std::sqrt(V[i]), 0.0,  r)  << std::endl;
                //stor.push_back(blackScholesCall(S[i], K,  T, std::sqrt(V[i]), 0.0,  r));
                stor.push_back(blackScholesPut(S[i], K,  T, std::sqrt(V[i]), 0.0,  r));
                vola.push_back(std::sqrt(V[i]));
                under.push_back(S[i]);
                T -= dt;
            }
        }
	    std::cout << std::accumulate(stor.begin() , stor.end(), 0.0)/stor.size() << std::endl;	
	    std::cout << std::accumulate(vola.begin() , vola.end(), 0.0)/vola.size() << std::endl;	
    }
    return 0;
}

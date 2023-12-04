#pragma once
// include everthing for this project here and always include this file
#include<iostream>
#include<algorithm>
#include<cmath>
#include<utility>
#include<memory>
#include<vector>
#include<numeric>
#include<fstream>
#include<iomanip>
#include<math.h>
using namespace std;

# define PI 3.14159265358979323846
# define E 2.71828182845904523536

enum Type { American, European };
enum CP { put, call };

double cdf(double d) {
	return (1 + std::erf(d / std::sqrt(2.0))) / 2.0;
}

class Option {
public:
	double S0;
	double K;
	double T;
	double sigma;
	double q;
	double r;
	Type type;
	CP cp;

	Option() {}
	Option(Type type_, CP cp_, double S0_, double K_, double T_, double sigma_, double q_,
		double r_) :S0(S0_), K(K_), T(T_), sigma(sigma_), q(q_), r(r_),
		type(type_), cp(cp_) {
		calculateBS();
	}

	double payoff(double S) {
		if (cp == put) {
			return max(K - S, 0.0);
		}
		else if (cp == call) {
			return max(S - K, 0.0);
		}
		return -1;
	}

	double BSvalue, BSdelta, BSvega;

	void calculateBS() {
		double d1 = (std::log(S0 / K) + (r - q + sigma * sigma / 2)*T) / sigma / std::sqrt(T);
		double d2 = d1 - sigma * std::sqrt(T);
		if (cp == call) {
			double Nd1 = cdf(d1);
			double Nd2 = cdf(d2);
			BSvalue = S0 * std::exp(-q * T)*Nd1 - K * std::exp(-r * T)*Nd2;
			BSdelta = std::exp(-q * T)*Nd1;
			//BSvega = S0 * exp(-q * T)*sqrt(T) / sqrt(2 * PI) * exp(-d1 * d1 / 2);
		}
		else if (cp == put) {
			double Nd1 = cdf(-d1);
			double Nd2 = cdf(-d2);
			BSvalue = -S0 * std::exp(-q * T)*Nd1 + K * std::exp(-r * T)*Nd2;
			BSdelta = -std::exp(-q * T)*Nd1;
			//BSvega = S0 * exp(-q * T)*sqrt(T) / sqrt(2 * PI) * exp(-d1 * d1 / 2);
		}	
	}


};



class DAO_Call_Euro_Option : public Option {
public:
    double B;
    
    DAO_Call_Euro_Option() {};
    DAO_Call_Euro_Option(Type European, CP call, double S0_, double K_, double T_, double sigma_, double q_, double r_, double B_) : B(B_){
        S0 = S0_;
        K = K_;
        T = T_;
        sigma = sigma_;
        q = q_;
        r = r_;
        type = European;
        cp = call;
        double s_temp = (1.0 * B * B) / S0;
        double a = (r-q)/sigma/sigma - 0.5;
        Option c1(European, call, S0, K, T, sigma, q, r);
        Option c2(European, call, s_temp, K, T, sigma, q, r);
        BSvalue = c1.BSvalue -  pow(((1.0 * B)/S0), 2.0 * a)* c2.BSvalue;
    }
    

};

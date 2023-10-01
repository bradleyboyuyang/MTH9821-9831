//
//  EuropeanOption.hpp
//  BinomialTree
//
//  Created by MS Wang on 2022/6/7.
//

#ifndef EuropeanOption_hpp
#define EuropeanOption_hpp

#include <functional>

class EuropeanOption {
public:
    double t_;      // Current time
    double S_;      // Spot price
    double K_;      // Strike price
    double T_;      // Maturity
    double sigma_;  // Volatility
    double r_;      // Const interest rate
    double q_;      // Dividend rate
    
private:
    // Intermediate values for Black-Scholes pricing
    double d1_;
    double d2_;
    double Zd1_;
    double Zd2_;
    double Nd1_;
    double Nd2_;
    double q_disc_; // Discount by dividend rate: exp(-q(T-t))
    double r_disc_; // Discount by interest rate: exp(-r(T-t))
    
    double Z(double t) const;   // PDF of a standard normal variable
    double Phi(double t) const; // CDF of a standard normal variable
    
    
public:
    EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q);
    ~EuropeanOption() = default;
    
    std::function<double (double, double)> CallPayoff() const;
    std::function<double (double, double)> PutPayoff() const;
    
    // Black-Scholes Price
    double Call() const;
    double Put() const;
    
    // Black-Scholes Greeks
    double DeltaCall() const;
    double DeltaPut() const;
    double GammaCall() const;
    double GammaPut() const;
    double ThetaCall() const;
    double ThetaPut() const;
    double VegaCall() const;
    double VegaPut() const;
    double RhoCall() const;
    double RhoPut() const;
};

#endif /* EuropeanOption_hpp */

//
//  TTree.hpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#ifndef TTree_hpp
#define TTree_hpp

#include <vector>
#include <functional>
#include <deque>
#include "EuropeanOption.hpp"
#include "TreePricer.hpp"

class TTree {
    // Trinomial tree pricer
    // Underlying asset follows a lognormal distribution
    // Continuous dividend
    // PI: path independent
    // EE: early exercise
protected:
    // Constructor takes in those values
    double S0_;     // Spot price at t = 0
    double sigma_;  // Volatility
    double T_;      // Maturity
    double steps_;  // Steps
    double r_;      // Interest rate
    double q_;      // Dividend rate
    
    // Helper values
    double dt_;     // delta t = T / steps
    double u_;      // u = exp(sigma * sqrt(2 * dt))
    double d_;      // d = 1 / u
    double r_disc_; // exp(-r*dt)
    double pu_;     // ((exp((r - q) * dt / 2) - exp(-sigma * sqrt(dt / 2))) / (exp(sigma * sqrt(dt / 2)) - exp(-sigma * sqrt(dt / 2)))) ^ 2
    double pd_;     // ((exp(sigma * sqrt(dt / 2)) - exp((r - q) * dt / 2)) / (exp(sigma * sqrt(dt / 2)) - exp(-sigma * sqrt(dt / 2)))) ^ 2
    double disc_pu; // disc * pu
    double disc_pd; // disc * pd
    double disc_pm; // disc * (1 - pu - pd)
    
public:
    TTree(double S0, double sigma, double T, std::size_t steps, double r, double q);
    TTree(const EuropeanOption& option, std::size_t steps);
    ~TTree() = default;
    
    TreeResult PathIndependentOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const;
    TreeResult EarlyExerciseOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const;
    
protected:
    std::deque<double> GenerateSMesh(double S0, double u, double d, std::size_t steps) const;
    
    std::vector<double> GeneratePayoff(const std::deque<double>& S, const std::function<double (double, double)>& payoff, double t) const;
    std::vector<double> GeneratePayoff(const std::deque<double> &S, const std::function<double (double, double)> &payoff, double t, const Barrier& barrier, double B) const;
    
    void EarlyExUpdate(std::vector<double>& V, const std::vector<double>& earlyex_payoff) const;
    
    void BacktrackPI(std::vector<double>& V_mesh) const;
    void BacktrackEE(std::vector<double>& V_mesh, std::deque<double>& S_mesh, size_t curr_step, const std::function<double (double, double)>& payoff, double curr_time) const;
    void BacktrackBarrier(std::vector<double>& V, std::deque<double>& S, const Barrier& barrier, double B) const;
    
    TreeResult GenTreeResult(double V00, double V10, double V11, double V12, double V20, double V22, double V24) const;
    
    TreeResult PIVanilla(const std::function<double (double, double)>& payoff) const;
//    TreeResult PIAvg(const std::function<double (double, double)>& payoff) const;
    TreeResult PIBS(const std::function<double (double, double)>& payoff) const;        // Black-Scholes
    TreeResult PIBSR(const std::function<double (double, double)>& payoff) const;       // Black-Scholes with Richardson extrapolation
    
    TreeResult EEVanilla(const std::function<double (double, double)>& payoff) const;
//    TreeResult EEAvg(const std::function<double (double, double)>& payoff) const;
    TreeResult EEBS(const std::function<double (double, double)>& payoff) const;        // Black-Scholes
    TreeResult EEBSR(const std::function<double (double, double)>& payoff) const;       // Black-Scholes with Richardson extrapolation
    
public:
    TreeResult BarrierVanilla(const std::function<double (double, double)>& payoff, const Barrier& barrier, double B) const;

};

#endif /* TTree_hpp */

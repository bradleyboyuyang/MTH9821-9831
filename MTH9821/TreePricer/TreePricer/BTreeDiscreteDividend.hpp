//
//  BTreeDiscreteDividend.hpp
//  TreePricer
//
//  Created by 王明森 on 11/23/22.
//

#ifndef BTreeDiscreteDividend_hpp
#define BTreeDiscreteDividend_hpp

#include "BTree.hpp"

class BTreeDiscDiv : public BTree {
    // Binomial tree pricer
    // Underlying asset follows a lognormal distribution
    // Allows discrete dividends
    // PI: path independent
    // EE: early exercise
protected:
    Dividend proportional_;
    Dividend fixed_;
    double equivalent_S0_;
    
public:
    BTreeDiscDiv(double S0, double sigma, double T, std::size_t steps, double r, double q, const Dividend& proportional, const Dividend& fixed);
    BTreeDiscDiv(const EuropeanOption& option, std::size_t steps, const Dividend& proportional, const Dividend& fixed);
    ~BTreeDiscDiv() = default;
    
    virtual TreeResult PathIndependentOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const override;
    virtual TreeResult EarlyExerciseOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const override;
    
    double EquivalentS0() const;
    
private:
    static double FindEquivalentS0(double S0, const Dividend& proportional, const Dividend& fixed, double r);
    
    template <typename container>
    std::vector<double> GeneratePayoff(const container& S, const std::function<double (double, double)>& payoff, double t, double cum_fixed) const ;
    
    TreeResult GenTreeResult(double V00, double V10, double V11, double V20, double V21, double V22, double S10, double S11, double S20, double S21, double S22) const;
    
    void BacktrackEE(std::vector<double>& V_mesh, std::array<std::deque<double>, 2>& S_mesh, size_t curr_step, const std::function<double (double, double)>& payoff, double curr_time, double cum_fixed) const;
    
    void ApplyDividend(std::array<std::deque<double>, 2>& S, size_t curr_step, double& cum_fixed, std::vector<std::size_t>& proportional_steps, std::vector<double>& proportional_dividends, std::vector<std::size_t>& fixed_steps, std::vector<double>& fixed_dividends) const;
    
    virtual TreeResult EEVanilla(const std::function<double (double, double)>& payoff) const override;
//    virtual TreeResult EEAvg(const std::function<double (double, double)>& payoff) const override;
//    virtual TreeResult EEBS(const std::function<double (double, double)>& payoff) const override;        // Black-Scholes
//    virtual TreeResult EEBSR(const std::function<double (double, double)>& payoff) const override;       // Black-Scholes with Richardson extrapolation
    
    std::vector<std::size_t> DateToStep(const std::vector<double>& dates, std::size_t steps, double T) const;
};

#endif /* BTreeDiscreteDividend_hpp */

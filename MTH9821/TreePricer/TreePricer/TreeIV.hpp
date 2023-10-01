//
//  TreeIV.hpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#ifndef TreeIV_hpp
#define TreeIV_hpp

#include "EuropeanOption.hpp"
#include "TreePricer.hpp"

class TreeIV {
private:
    double S_;      // Spot price
    double K_;      // Strike price
    double T_;      // Maturity
    double r_;      // Const interest rate
    double q_;      // Dividend rate
    double V0_;     // Value at t=0
    
    static double tol_step;     // Tolerance used in optimal step determination
    static double tol_approx;   // Tolerance used in approximation
    static std::size_t step_init;   // Optimal step initiation
    static double vol_init;     
    
public:
    TreeIV(double S, double K, double T, double r, double q, double V0);
    TreeIV(const EuropeanOption& option);
    ~TreeIV() = default;
    
    double BTreeEuropeanPut(const TreeModifier& modifier) const;
    double BTreeAmericanPut(const TreeModifier& modifier) const;
    double TTreeEuropeanPut(const TreeModifier& modifier) const;
    double TTreeAmericanPut(const TreeModifier& modifier) const;
    
private:
    std::size_t BTreePutSteps(const TreeModifier& modifier) const;
    std::size_t TTreePutSteps(const TreeModifier& modifier) const;
};

#endif /* TreeIV_hpp */

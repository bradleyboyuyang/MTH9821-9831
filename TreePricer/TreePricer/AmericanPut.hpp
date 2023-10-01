//
//  AmericanPut.hpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#ifndef AmericanPut_hpp
#define AmericanPut_hpp

#include "TreePricer.hpp"
#include "EuropeanOption.hpp"

class AmericanPut {
private:
    EuropeanOption option_;
    
public:
    AmericanPut(double S, double K, double T, double sigma, double r, double q);
    AmericanPut(const EuropeanOption& option);
    
    // One-time pricer
    TreeResult BinomialTree(std::size_t steps, const TreeModifier& modifier, bool variance_reduction = false) const;
    TreeResult TrinomialTree(std::size_t steps, const TreeModifier& modifier, bool variance_reduction = false) const;
    TreeResult BinomialTree(std::size_t steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed, bool variance_reduction = false) const;
    
    // Iterative pricer
    std::tuple<TreeResult, std::size_t> BinomialTree(std::size_t init_steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed, double tol, bool variance_reduction = false) const;
};

#endif /* AmericanPut_hpp */

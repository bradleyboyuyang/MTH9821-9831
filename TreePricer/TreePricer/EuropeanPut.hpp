//
//  EuropeanPut.hpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#ifndef EuropeanPut_hpp
#define EuropeanPut_hpp

#include "TreePricer.hpp"
#include "EuropeanOption.hpp"
#include <tuple>

class EuropeanPut {
private:
    EuropeanOption option_;
    
public:
    EuropeanPut(double S, double K, double T, double sigma, double r, double q);
    EuropeanPut(const EuropeanOption& option);
    
    // One-time pricer
    TreeResult BinomialTree(std::size_t steps, const TreeModifier& modifier) const;
    TreeResult TrinomialTree(std::size_t steps, const TreeModifier& modifier) const;
    TreeResult BinomialTree(std::size_t steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed) const;
    
    // Iterative pricer
    std::tuple<TreeResult, std::size_t> BinomialTree(std::size_t init_steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed, double tol) const;
    
    TreeResult TrinomialTree(std::size_t steps, const TreeModifier& modifier, const Barrier& barrier, double B) const;
    
    TreeResult BinomialTree(std::size_t steps, const TreeModifier& modifier, const Barrier& barrier, double B) const;
};

#endif /* EuropeanPut_hpp */

//
//  EuropeanPut.cpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#include "EuropeanPut.hpp"
#include "BTree.hpp"
#include "TTree.hpp"
#include "BTreeDiscreteDividend.hpp"

EuropeanPut::EuropeanPut(double S, double K, double T, double sigma, double r, double q) : option_(0., S, K, T, sigma, r, q) {}
EuropeanPut::EuropeanPut(const EuropeanOption& option) : option_(option) {}

TreeResult EuropeanPut::BinomialTree(std::size_t steps, const TreeModifier &modifier) const {
    
    BTree european_tree(option_, steps);
    
    TreeResult res = european_tree.PathIndependentOption(option_.PutPayoff(), modifier);
    
    return res;
}

TreeResult EuropeanPut::TrinomialTree(std::size_t steps, const TreeModifier &modifier) const {
    
    TTree european_tree(option_, steps);
    
    TreeResult res = european_tree.PathIndependentOption(option_.PutPayoff(), modifier);
    
    return res;
}

TreeResult EuropeanPut::BinomialTree(std::size_t steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed) const {
    
    BTreeDiscDiv european_tree(option_, steps, proportional, fixed);
    
    TreeResult res = european_tree.PathIndependentOption(option_.PutPayoff(), modifier);
    
    return res;
}

// Iterative pricer
std::tuple<TreeResult, std::size_t> EuropeanPut::BinomialTree(std::size_t init_steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed, double tol) const {
    std::size_t steps = init_steps;
    
    TreeResult res = this->BinomialTree(steps, modifier, proportional, fixed);
    double v_new = res.value;
    double v_old = -1. - tol;
    
    while (std::abs(v_new - v_old) >= tol) {
        steps <<= 1;
        v_old = v_new;
        res = this->BinomialTree(steps, modifier, proportional, fixed);
        v_new = res.value;
    }
    
    return std::make_tuple(res, steps);
}

TreeResult EuropeanPut::TrinomialTree(std::size_t steps, const TreeModifier& modifier, const Barrier& barrier, double B) const {
    
    // !!!: This is a DaO CALL!!!
    
    TTree european_tree(option_, steps);
    TreeResult res =  european_tree.BarrierVanilla(option_.CallPayoff(), DownAndOut, B);
    
    return res;
}

TreeResult EuropeanPut::BinomialTree(std::size_t steps, const TreeModifier& modifier, const Barrier& barrier, double B) const {
    BTree european_tree(option_, steps);
    TreeResult res = european_tree.BarrierVanilla(option_.CallPayoff(), DownAndOut, B);
    
    return res;
}

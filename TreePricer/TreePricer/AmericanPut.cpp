//
//  AmericanPut.cpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#include "AmericanPut.hpp"
#include "BTree.hpp"
#include "TTree.hpp"
#include "BTreeDiscreteDividend.hpp"

AmericanPut::AmericanPut(double S, double K, double T, double sigma, double r, double q) : option_(0., S, K, T, sigma, r, q) {}
AmericanPut::AmericanPut(const EuropeanOption& option) : option_(option) {}

TreeResult AmericanPut::BinomialTree(std::size_t steps, const TreeModifier &modifier, bool variance_reduction) const {
    
    BTree tree(option_, steps);
    
    TreeResult res = tree.EarlyExerciseOption(option_.PutPayoff(), modifier);
    
    if (variance_reduction) {
        TreeResult euro_res = tree.PathIndependentOption(option_.PutPayoff(), modifier);
        
        res.value -= (euro_res.value - option_.Put());
        res.delta -= (euro_res.delta - option_.DeltaPut());
        res.gamma -= (euro_res.gamma - option_.GammaPut());
        res.theta -= (euro_res.theta - option_.ThetaPut());
    }
    
    return res;
}

TreeResult AmericanPut::TrinomialTree(std::size_t steps, const TreeModifier &modifier, bool variance_reduction) const {
    
    TTree tree(option_, steps);
    
    TreeResult res = tree.EarlyExerciseOption(option_.PutPayoff(), modifier);
    
    if (variance_reduction) {
        TreeResult euro_res = tree.PathIndependentOption(option_.PutPayoff(), modifier);
        
        res.value -= (euro_res.value - option_.Put());
        res.delta -= (euro_res.delta - option_.DeltaPut());
        res.gamma -= (euro_res.gamma - option_.GammaPut());
        res.theta -= (euro_res.theta - option_.ThetaPut());
    }
    
    return res;
}

TreeResult AmericanPut::BinomialTree(std::size_t steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed, bool variance_reduction) const {
    
    BTreeDiscDiv tree(option_, steps, proportional, fixed);
    
    TreeResult res = tree.EarlyExerciseOption(option_.PutPayoff(), modifier);
    
    if (variance_reduction) {
        TreeResult euro_res = tree.PathIndependentOption(option_.PutPayoff(), modifier);
        EuropeanOption equivalent_option(0., tree.EquivalentS0(), option_.K_, option_.T_, option_.sigma_, option_.r_, option_.q_);
        
        res.value -= (euro_res.value - equivalent_option.Put());
        res.delta -= (euro_res.delta - equivalent_option.DeltaPut());
        res.gamma -= (euro_res.gamma - equivalent_option.GammaPut());
        res.theta -= (euro_res.theta - equivalent_option.ThetaPut());
    }
    
    return res;
    
}

// Iterative pricer
std::tuple<TreeResult, std::size_t> AmericanPut::BinomialTree(std::size_t init_steps, const TreeModifier& modifier, const Dividend& proportional, const Dividend& fixed, double tol, bool variance_reduction) const {
    std::size_t steps = init_steps;
    
    TreeResult res = this->BinomialTree(steps, modifier, proportional, fixed, variance_reduction);
    double v_new = res.value;
    double v_old = -1. - tol;
    
    while (std::abs(v_new - v_old) >= tol) {
        steps <<= 1;
        v_old = v_new;
        res = this->BinomialTree(steps, modifier, proportional, fixed, variance_reduction);
        v_new = res.value;
    }
    
    return std::make_tuple(res, steps);
}

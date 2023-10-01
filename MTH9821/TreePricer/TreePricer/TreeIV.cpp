//
//  TreeIV.cpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#include "TreeIV.hpp"
#include "BTree.hpp"
#include "TTree.hpp"
#include <cmath>
#include "EuropeanPut.hpp"
#include "AmericanPut.hpp"

double TreeIV::tol_step = std::pow(10, -6);     // Tolerance used in optimal step determination
double TreeIV::tol_approx = std::pow(10, -6);   // Tolerance used in approximation
std::size_t TreeIV::step_init = 10;             // Optimal step initiation
double TreeIV::vol_init = .25;

TreeIV::TreeIV(double S, double K, double T, double r, double q, double V0) : S_(S), K_(K), T_(T), r_(r), q_(q), V0_(V0) {}

std::size_t TreeIV::BTreePutSteps(const TreeModifier& modifier) const {
    
    std::size_t steps = step_init >> 1;
    
    EuropeanOption option(0., S_, K_, T_, vol_init, r_, q_);
    EuropeanPut put(option);
    
    double tree_estimate = option.Put() + 1.;
    double tree_estimate_old = 0.;
    
    while (std::abs(tree_estimate - tree_estimate_old) >= tol_step) {
        steps <<= 1;
        tree_estimate_old = tree_estimate;
        tree_estimate = put.BinomialTree(steps, modifier).value;
        std::cout << steps << std::endl;
    }
    
    return steps;
}

std::size_t TreeIV::TTreePutSteps(const TreeModifier& modifier) const {
    
    std::size_t steps = step_init >> 1;
    
    EuropeanOption option(0., S_, K_, T_, vol_init, r_, q_);
    EuropeanPut put(option);
    
    double tree_estimate = option.Put() + 1.;
    double tree_estimate_old = 0.;
    
    while (std::abs(tree_estimate - tree_estimate_old) >= tol_step) {
        steps <<= 1;
        tree_estimate_old = tree_estimate;
        tree_estimate = put.TrinomialTree(steps, modifier).value;
        std::cout << steps << std::endl;
    }
    
    return steps;
}

double TreeIV::BTreeEuropeanPut(const TreeModifier& modifier) const {
    
    std::size_t steps = this->BTreePutSteps(modifier);
    
    double sig_last = .05;
    double sig_new = 1.;
    
    while (std::abs(sig_new - sig_last) >= tol_approx) {
        EuropeanPut new_put(S_, K_, T_, sig_new, r_, q_);
        EuropeanPut last_put(S_, K_, T_, sig_last, r_, q_);
        
        double v_new = new_put.BinomialTree(steps, modifier).value;
        double v_last = last_put.BinomialTree(steps, modifier).value;
        
        // Secant method
        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);
        
        sig_last = sig_new;
        sig_new = sig_newest;
    }
    
    return sig_new;
}

double TreeIV::TTreeEuropeanPut(const TreeModifier& modifier) const {
    
    std::size_t steps = this->TTreePutSteps(modifier);
    
    double sig_last = .05;
    double sig_new = 1.;
    
    while (std::abs(sig_new - sig_last) >= tol_approx) {
        EuropeanPut new_put(S_, K_, T_, sig_new, r_, q_);
        EuropeanPut last_put(S_, K_, T_, sig_last, r_, q_);
        
        double v_new = new_put.TrinomialTree(steps, modifier).value;
        double v_last = last_put.TrinomialTree(steps, modifier).value;
        
        // Secant method
        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);
        
        sig_last = sig_new;
        sig_new = sig_newest;
    }
    
    return sig_new;
}

double TreeIV::BTreeAmericanPut(const TreeModifier& modifier) const {
    std::size_t steps = this->BTreePutSteps(modifier);
    
    double sig_last = .05;
    double sig_new = 1.;
    
    while (std::abs(sig_new - sig_last) >= tol_approx) {
        AmericanPut new_put(S_, K_, T_, sig_new, r_, q_);
        AmericanPut last_put(S_, K_, T_, sig_last, r_, q_);
        
        double v_new = new_put.BinomialTree(steps, modifier).value;
        double v_last = last_put.BinomialTree(steps, modifier).value;
        
        // Secant method
        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);
        
        sig_last = sig_new;
        sig_new = sig_newest;
    }
    
    return sig_new;
}

double TreeIV::TTreeAmericanPut(const TreeModifier& modifier) const {
    std::size_t steps = this->TTreePutSteps(modifier);
    
    double sig_last = .05;
    double sig_new = 1.;
    
    while (std::abs(sig_new - sig_last) >= tol_approx) {
        AmericanPut new_put(S_, K_, T_, sig_new, r_, q_);
        AmericanPut last_put(S_, K_, T_, sig_last, r_, q_);
        
        double v_new = new_put.TrinomialTree(steps, modifier).value;
        double v_last = last_put.TrinomialTree(steps, modifier).value;
        
        // Secant method
        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);
        
        sig_last = sig_new;
        sig_new = sig_newest;
    }
    
    return sig_new;
}

//
//  BTreeDiscreteDividend.cpp
//  TreePricer
//
//  Created by 王明森 on 11/23/22.
//

#include "BTreeDiscreteDividend.hpp"
#include <cmath>
#include<cassert>

BTreeDiscDiv::BTreeDiscDiv(double S0, double sigma, double T, std::size_t steps, double r, double q, const Dividend& proportional, const Dividend& fixed) : BTree(S0, sigma, T, steps, r, q), proportional_(proportional), fixed_(fixed) {
    assert(q == 0.);
    equivalent_S0_ = BTreeDiscDiv::FindEquivalentS0(S0, proportional, fixed, r);
}

BTreeDiscDiv::BTreeDiscDiv(const EuropeanOption& option, std::size_t steps, const Dividend& proportional, const Dividend& fixed) : BTree(option, steps), proportional_(proportional), fixed_(fixed) {
    assert(option.q_ == 0.);
    equivalent_S0_ = BTreeDiscDiv::FindEquivalentS0(option.S_, proportional, fixed, option.r_);
}

double BTreeDiscDiv::FindEquivalentS0(double S0, const Dividend& proportional, const Dividend& fixed, double r) {
    
    double S = S0;
    // The "sophisticated" algorithm
    auto tit_prop = proportional.dates.cbegin();
    auto dit_prop = proportional.dividends.cbegin();
    auto tit_fixed = fixed.dates.cbegin();
    auto dit_fixed = fixed.dividends.cbegin();

    while ((tit_prop != proportional.dates.cend()) || (tit_fixed != fixed.dates.cend())) {
        if (tit_prop == proportional.dates.cend()) {
            // Proportional dividends exhausted, apply all remaining fixed dividends
            while (tit_fixed != fixed.dates.cend()) {
                S -= (*dit_fixed) * std::exp(-r * (*tit_fixed));
                dit_fixed++;
                tit_fixed++;
            }
            break;
        } else if (tit_fixed == fixed.dates.cend()) {
            // Fixed dividends exhausted, apply all remaining proportional dividends
            while (tit_prop != proportional.dates.cend()) {
                S *= (1. - *dit_prop);
                dit_prop++;
                tit_prop++;
            }
            break;
        } else {
            // Both types of dividends expected. Apply the earliest dividend possible.
            // Proportional dividends go first if (in the unlikely case that) two dividends coincide in time
            if ((*tit_fixed) < (*tit_prop)) {
                // The fixed dividend is earlier
                S -= (*dit_fixed) * std::exp(-r * (*tit_fixed));
                dit_fixed++;
                tit_fixed++;
            } else {
                // The proportional dividend is earlier
                S *= (1. - *dit_prop);
                dit_prop++;
                tit_prop++;
            }
        }
    }

//    // The "naive" algorithm
//    for (auto dit = proportional.dividends.cbegin(); dit != proportional.dividends.cend(); dit++) {
//        S *= (1. - *dit);
//    }
//
//    for (auto tit = fixed.dates.cbegin(), dit = fixed.dividends.cbegin(); tit != fixed.dates.cend() && dit != fixed.dividends.cend(); tit++, dit++) {
//        S -= (*dit) * std::exp(-r * (*tit));
//    }
    
    return S;
}

double BTreeDiscDiv::EquivalentS0() const {
    return equivalent_S0_;
}

TreeResult BTreeDiscDiv::PathIndependentOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const {
    // Since the option is path independent, we can simply construct a continuous-dividend binomial tree with adjusted S0 for pricing.
    BTree continuous_dividend_tree(equivalent_S0_, sigma_, T_, steps_, r_, q_);
    
    // Adjust Greeks
    double multiplier = 1.;
    for (double div : proportional_.dividends) {
        multiplier *= (1. - div);
    }
    TreeResult continuous_res(continuous_dividend_tree.PathIndependentOption(payoff, modifier));
    
    continuous_res.delta *= multiplier;
    continuous_res.gamma *= multiplier * multiplier;
    
    return continuous_res;
}

TreeResult BTreeDiscDiv::EarlyExerciseOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const {
    // We need to rewrite almost everything since the tree itself is perturbed.
    // Only the vanilla version is needed.
    switch (modifier) {
        case TreeModifier::vanilla:
            // Vanilla tree
            return this->EEVanilla(payoff);
            break;
//        case TreeModifier::average:
//            // Average tree
//            return this->EEAvg(payoff);
//            break;
//        case TreeModifier::BS:
//            // Black-Scholes tree
//            return this->EEBS(payoff);
//            break;
//        case TreeModifier::BSR:
//            // Black-Scholes tree with Richardson extrapolation
//            return this->EEBSR(payoff);
//            break;
        default:
            // Default to vanilla tree
            return this->EEVanilla(payoff);
            break;
    }
}

std::vector<std::size_t> BTreeDiscDiv::DateToStep(const std::vector<double>& dates, std::size_t steps, double T) const {
    auto converter = [=](double date)->std::size_t {
        return static_cast<std::size_t>(std::round(date * steps / T));
    };
    
    std::vector<std::size_t> res(dates.size());
    std::transform(dates.cbegin(), dates.cend(), res.begin(), converter);
    return res;
}

template <typename container>
std::vector<double> BTreeDiscDiv::GeneratePayoff(const container& S, const std::function<double (double, double)>& payoff, double t, double cum_fixed) const {
    auto payoff_t = std::bind(payoff, std::placeholders::_1, t);
    
//    for (double s : S) {
//        std::cout << s << '\t';
//    }
//    std::cout << std::endl;
    
    std::vector<double> V(S.cbegin(), S.cend());
    std::for_each(V.begin(), V.end(), [=](double& s)->void { s += cum_fixed; });
    std::transform(V.cbegin(), V.cend(), V.begin(), payoff_t);
    return V;
}

TreeResult BTreeDiscDiv::GenTreeResult(double V00, double V10, double V11, double V20, double V21, double V22, double S10, double S11, double S20, double S21, double S22) const {
    double delta = (V10 - V11) / (S10 - S11);
    double gamma = ((V20 - V21) / (S20 - S21) - (V21 - V22) / (S21 - S22)) / (.5 * (S20 - S22));
    double theta = (V21 - V00) / (2. * dt_);
    
    return TreeResult({V00, delta, gamma, theta});
}

void BTreeDiscDiv::BacktrackEE(std::vector<double>& V_mesh, std::array<std::deque<double>, 2>& S_mesh, size_t curr_step, const std::function<double (double, double)>& payoff, double curr_time, double cum_fixed) const {
    
    // Backtrack as if path independent
    this->BTree::BacktrackPI(V_mesh);
    
    // We maintain two vectors of S and we use them alternatingly.
    // The current S[i] always has the same size as V.
    // The other has (curr_step) elements.
    size_t which_S = (S_mesh[1].size() == (V_mesh.size()));
    
    // Find the payoff if we exercise the option at this step
    std::vector<double> earlyex_payoff(this->GeneratePayoff(S_mesh[which_S], payoff, curr_time, cum_fixed));
    
    // If early exercise is more profitable, replace value with early exercise payoff.
    this->BTree::EarlyExUpdate(V_mesh, earlyex_payoff);
    
    // Also shrink the current S since we will have less steps in the next two iterations.
    if (S_mesh[which_S].size() >= 2) {
        S_mesh[which_S].pop_back();
        S_mesh[which_S].pop_front();
    }
}

void BTreeDiscDiv::ApplyDividend(std::array<std::deque<double>, 2>& S, size_t curr_step, double& cum_fixed, std::vector<std::size_t>& proportional_steps, std::vector<double>& proportional_dividends, std::vector<std::size_t>& fixed_steps, std::vector<double>& fixed_dividends) const {
    
    // Discount existing cumulative fixed dividend.
    cum_fixed *= std::exp(-r_ * dt_);
    
    // Check fixed dividend
    while ((!fixed_steps.empty()) && (fixed_steps.back() == curr_step)) {
        // Apply dividend
        cum_fixed += fixed_dividends.back();
        
        // Shorten fixed dividend array
        fixed_steps.pop_back();
        fixed_dividends.pop_back();
    }
    
    // Check proportional dividend
    while ((!proportional_steps.empty()) && (proportional_steps.back() == curr_step)) {
        // Apply dividend
        std::for_each(S[0].begin(), S[0].end(), [=](double& s)->void { s /= (1. - proportional_dividends.back()); });
        std::for_each(S[1].begin(), S[1].end(), [=](double& s)->void { s /= (1. - proportional_dividends.back()); });
        cum_fixed /= (1. - proportional_dividends.back());
        // Shorten proportional dividend array
        proportional_steps.pop_back();
        proportional_dividends.pop_back();
    }
}

TreeResult BTreeDiscDiv::EEVanilla(const std::function<double (double, double)>& payoff) const {
    // Convert dividend dates to steps (to avoid rounding errors)
    std::vector<std::size_t> proportional_steps(this->DateToStep(proportional_.dates, steps_, T_));
    std::vector<std::size_t> fixed_steps(this->DateToStep(fixed_.dates, steps_, T_));
    
    // Copy the dividend amounts so that we can pop from the back.
    std::vector<double> proportional_dividends(proportional_.dividends);
    std::vector<double> fixed_dividends(fixed_.dividends);
    
    // Current cumulative fixed dividend (starting from the end)
    double cum_fixed = 0.;
    
    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::array<std::deque<double>, 2> S(this->BTree::GenerateSMeshEE(equivalent_S0_, u_, d_, steps_));
    
    double curr_time = T_;
    std::size_t curr_step = steps_;
    std::vector<double> V(this->GeneratePayoff(S[0], payoff, curr_time, cum_fixed));
    S[0].pop_back();
    S[0].pop_front();
    this->ApplyDividend(S, curr_step, cum_fixed, proportional_steps, proportional_dividends, fixed_steps, fixed_dividends);
    
    curr_step--;
    curr_time -= dt_;
    
    double S20 = 0., S21 = 0., S22 = 0.;
    std::size_t which_S;
    while (curr_step >= 2) {
        // During backtracking, curr_step and curr_time are exactly the current step (0-indexed) and the current time
        if (curr_step == 2) {
            which_S = std::size_t(S[1].size() == 3);
            assert(S[which_S].size() == 3);
            S20 = S[which_S][0] + cum_fixed;
            S21 = S[which_S][1] + cum_fixed;
            S22 = S[which_S][2] + cum_fixed;
        }
        this->BacktrackEE(V, S, curr_step, payoff, curr_time, cum_fixed);
        this->ApplyDividend(S, curr_step, cum_fixed, proportional_steps, proportional_dividends, fixed_steps, fixed_dividends);
        curr_step--;
        curr_time -= dt_;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];
    assert(V.size() == 3);
    
    // Backtrack
    which_S = std::size_t(S[1].size() == 2);
    assert(S[which_S].size() == 2);
    double S10 = S[which_S][0] + cum_fixed;
    double S11 = S[which_S][1] + cum_fixed;
    this->BacktrackEE(V, S, curr_step, payoff, curr_time, cum_fixed);
    this->ApplyDividend(S, curr_step, cum_fixed, proportional_steps, proportional_dividends, fixed_steps, fixed_dividends);
    curr_step--;
    curr_time -= dt_;
    
    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time, cum_fixed);
    this->ApplyDividend(S, curr_step, cum_fixed, proportional_steps, proportional_dividends, fixed_steps, fixed_dividends);
    curr_step--;
    curr_time -= dt_;
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V20, V21, V22, S10, S11, S20, S21, S22);
}

#ifndef BTree_cpp
#define BTree_cpp

#include "BTree.hpp"
#include <cmath>
#include <algorithm>
#include<cassert>
#include <iostream>

BTree::BTree(double S0, double sigma, double T, std::size_t steps, double r, double q) : S0_(S0), sigma_(sigma), T_(T), steps_(steps), r_(r), q_(q),
    dt_(T / steps),
    u_(std::exp(sigma * std::sqrt(dt_))),
    d_(1. / u_),
    r_disc_(std::exp(-r * dt_)),
    p_((std::exp((r - q) * dt_) - d_) / (u_ - d_)),
    disc_p(r_disc_ * p_),
    disc_1p(r_disc_ * (1. - p_))
{}

BTree::BTree(const EuropeanOption& option, std::size_t steps) : BTree(option.S_, option.sigma_, option.T_, steps, option.r_, option.q_) {}

TreeResult BTree::PathIndependentOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const {
    switch (modifier) {
        case TreeModifier::vanilla:
            // Vanilla tree
            return this->PIVanilla(payoff);
            break;
        case TreeModifier::average:
            // Average tree
            return this->PIAvg(payoff);
            break;
        case TreeModifier::BS:
            // Black-Scholes tree
            return this->PIBS(payoff);
            break;
        case TreeModifier::BSR:
            // Black-Scholes tree with Richardson extrapolation
            return this->PIBSR(payoff);
            break;
        default:
            // Default to vanilla tree
            return this->PIVanilla(payoff);
            break;
    }
}

TreeResult BTree::EarlyExerciseOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const {
    switch (modifier) {
        case TreeModifier::vanilla:
            // Vanilla tree
            return this->EEVanilla(payoff);
            break;
        case TreeModifier::average:
            // Average tree
            return this->EEAvg(payoff);
            break;
        case TreeModifier::BS:
            // Black-Scholes tree
            return this->EEBS(payoff);
            break;
        case TreeModifier::BSR:
            // Black-Scholes tree with Richardson extrapolation
            return this->EEBSR(payoff);
            break;
        default:
            // Default to vanilla tree
            return this->EEVanilla(payoff);
            break;
    }
}

std::vector<double> BTree::GenerateSMeshPI(double S0, double u, double d, std::size_t steps) const {
    std::vector<double> S({S0 * std::pow(u, static_cast<double>(steps))});
    double d2 = d * d;
    for (std::size_t i = 0; i < steps; i++) {
        S.push_back(S.back() * d2);
    }
    return S;
}

std::array<std::deque<double>, 2> BTree::GenerateSMeshEE(double S0, double u, double d, std::size_t steps) const {
    
    std::array<std::deque<double>, 2> S;
    std::deque<double> longer_mesh;
    std::deque<double> shorter_mesh;
    
    longer_mesh.push_back(S0 * std::pow(u, static_cast<double>(steps)));
    shorter_mesh.push_back(S0 * std::pow(u, static_cast<double>(steps) - 1.));
    
    double d2 = d * d;
    for (std::size_t i = 0; i < steps - 1; i++) {
        longer_mesh.push_back(longer_mesh.back() * d2);
        shorter_mesh.push_back(shorter_mesh.back() * d2);
    }
    longer_mesh.push_back(longer_mesh.back() * d2);
    // Now, longer_mesh has steps + 1 elements and shorter_mesh has steps elements
    // The current S[i] always has (curr_step + 1) elements.

    return std::array<std::deque<double>, 2>({std::move(longer_mesh), std::move(shorter_mesh)});
}

template <typename container>
std::vector<double> BTree::GeneratePayoff(const container& S, const std::function<double (double, double)>& payoff, double t) const {
    auto payoff_t = std::bind(payoff, std::placeholders::_1, t);
    
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), payoff_t);
    
    return V;
}

std::vector<double> BTree::GeneratePayoffBarrier(const std::deque<double>& S, const std::function<double (double, double)>& payoff, double t, const Barrier& barrier, double B) const {
    // MARK: Barrier payoff
    
    assert(barrier == DownAndOut);
    
    auto payoff_t = std::bind(payoff, std::placeholders::_1, t);
    
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), payoff_t);
    
    std::transform(S.cbegin(), S.cend(), V.cbegin(), V.begin(), [&](double s, double v)->double {
        if (s <= B) {
            return 0.;
        } else {
            return v;
        }
    });
    
    return V;
}

void BTree::EarlyExUpdate(std::vector<double>& V, const std::vector<double>& earlyex_payoff) const {
    assert(V.size() == earlyex_payoff.size());
    
    auto earlyex_it = earlyex_payoff.cbegin();
    for (auto Vit = V.begin(); Vit != V.end(); Vit++) {
        *Vit = std::max(*Vit, *earlyex_it);
        earlyex_it++;
    }
}

void BTree::BacktrackPI(std::vector<double>& V_mesh) const {
    for (auto Vit = V_mesh.begin(); Vit + 1 != V_mesh.end(); Vit++) {
        // Get value of the last node by taking the expectation and discounting
        // V = disc(p_u * V_u + p_d * V_d)
        *Vit = disc_p * (*Vit) + disc_1p * (*(Vit + 1));
    }
    
    // Shrink V since there are fewer nodes needed.
    V_mesh.pop_back();
}

void BTree::BacktrackEE(std::vector<double>& V_mesh, std::array<std::deque<double>, 2>& S_mesh, size_t curr_step, const std::function<double (double, double)>& payoff, double curr_time) const {
    
    // Backtrack as if path independent
    this->BacktrackPI(V_mesh);
    
    // We maintain two vectors of S and we use them alternatingly.
    // The current S[i] always has (curr_step + 1) elements.
    // The other has (curr_step) elements.
    size_t which_S = (S_mesh[1].size() == (curr_step + 1));
    
    // Find the payoff if we exercise the option at this step
    std::vector<double> earlyex_payoff(this->GeneratePayoff(S_mesh[which_S], payoff, curr_time));
    
    // If early exercise is more profitable, replace value with early exercise payoff.
    this->EarlyExUpdate(V_mesh, earlyex_payoff);
    
    // Also shrink the current S since we will have less steps in the next two iterations.
    if (S_mesh[which_S].size() >= 2) {
        S_mesh[which_S].pop_back();
        S_mesh[which_S].pop_front();
    }
}

void BTree::BacktrackBarrier(std::vector<double>& V, std::array<std::deque<double>, 2>& S, size_t curr_step, const Barrier& barrier, double B) const {
    
    assert(barrier == DownAndOut);
    // Backtrack as if path independent
    this->BacktrackPI(V);
    
    // We maintain two vectors of S and we use them alternatingly.
    // The current S[i] always has (curr_step + 1) elements.
    // The other has (curr_step) elements.
    size_t which_S = (S[1].size() == (curr_step + 1));
    
    // Set to zero if S is down
    std::transform(S[which_S].cbegin(), S[which_S].cend(), V.cbegin(), V.begin(), [&](double s, double v)->double {
        if (s <= B) {
            return 0.;
        } else {
            return v;
        }
    });
    
    // Also shrink the current S since we will have less steps in the next two iterations.
    if (S[which_S].size() >= 2) {
        S[which_S].pop_back();
        S[which_S].pop_front();
    }
}

TreeResult BTree::GenTreeResult(double V00, double V10, double V11, double V20, double V21, double V22) const {
    
    double delta = (V10 - V11) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V21) / (S0_ * (u_ * u_ - 1.)) - ((V21 - V22) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V21 - V00)/ (2. * dt_);
    
    return TreeResult({V00, delta, gamma, theta});
}

TreeResult BTree::PIVanilla(const std::function<double (double, double)>& payoff) const {
    
    // Generate asset price at maturity
    std::vector<double> S(this->GenerateSMeshPI(S0_, u_, d_, steps_));
    
    // Generate derivative value at maturity
    std::vector<double> V(this->GeneratePayoff(S, payoff, T_));
    
    // Backtrack until t = 2 * dt (the resulting sub-tree is used for Greek calculation)
    std::size_t curr_step = steps_ - 1;
    while (curr_step >= 2) {
        this->BacktrackPI(V);
        curr_step--;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];
    assert(V.size() == 3);
    
    // Backtrack
    this->BacktrackPI(V);
    
    // Store nodes of Step 2
    double V10 = V[0];
    double V11 = V[1];
    
    // Backtrack
    this->BacktrackPI(V);
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V20, V21, V22);
}

TreeResult BTree::EEVanilla(const std::function<double (double, double)>& payoff) const {
    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::array<std::deque<double>, 2> S(this->GenerateSMeshEE(S0_, u_, d_, steps_));
    
    double curr_time = T_;
    std::size_t curr_step = steps_;
    std::vector<double> V(this->GeneratePayoff(S[0], payoff, curr_time));
    S[0].pop_back();
    S[0].pop_front();
    curr_step--;
    curr_time -= dt_;
    
    while (curr_step >= 2) {
        // During backtracking, curr_step and curr_time are exactly the current step (0-indexed) and the current time
        this->BacktrackEE(V, S, curr_step, payoff, curr_time);
        curr_step--;
        curr_time -= dt_;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];
    assert(V.size() == 3);
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_step--;
    curr_time -= dt_;
    
    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_step--;
    curr_time -= dt_;
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V20, V21, V22);
}

TreeResult BTree::PIAvg(const std::function<double (double, double)>& payoff) const {
    TreeResult res1 = this->PIVanilla(payoff);
    
    BTree longer_tree(S0_, sigma_, T_, steps_ + 1, r_, q_);
    
    TreeResult res2 = longer_tree.PIVanilla(payoff);
    
    return TreeResult({(res1.value + res2.value) / 2., (res1.delta + res2.delta) / 2., (res1.gamma + res2.gamma) / 2., (res1.theta + res2.theta) / 2.});
}

TreeResult BTree::EEAvg(const std::function<double (double, double)>& payoff) const {
    TreeResult res1 = this->EEVanilla(payoff);
    
    BTree longer_tree(S0_, sigma_, T_, steps_ + 1, r_, q_);
    
    TreeResult res2 = longer_tree.EEVanilla(payoff);
    
    return TreeResult({(res1.value + res2.value) / 2., (res1.delta + res2.delta) / 2., (res1.gamma + res2.gamma) / 2., (res1.theta + res2.theta) / 2.});
}

TreeResult BTree::PIBS(const std::function<double (double, double)>& payoff) const {
    std::size_t curr_step = steps_;
    curr_step--;    // We do not need the last step here
    
    // Generate asset price one step before maturity
    std::vector<double> S(this->GenerateSMeshPI(S0_, u_, d_, curr_step));
    
    // Generate the BS pricer
    auto last_step = [&](double S_T, double t)->double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);
        
        return option.Put();
    };
    
    
    // Generate derivative value at the step just before maturity
    std::vector<double> V(this->GeneratePayoff(S, last_step, 0.));  // The last param is a dummy
    curr_step--;
    
    // Backtrack until t = 2 * dt (the resulting sub-tree is used for Greek calculation)
    while (curr_step >= 2) {
        this->BacktrackPI(V);
        curr_step--;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];
    assert(V.size() == 3);
    
    // Backtrack
    this->BacktrackPI(V);
    
    // Store nodes of Step 2
    double V10 = V[0];
    double V11 = V[1];
    
    // Backtrack
    this->BacktrackPI(V);
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V20, V21, V22);
}

TreeResult BTree::EEBS(const std::function<double (double, double)>& payoff) const {
    std::size_t curr_step = steps_;
    double curr_time = T_;
    curr_step--;    // We do not need the last step here
    curr_time -= dt_;
    
    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::array<std::deque<double>, 2> S(this->GenerateSMeshEE(S0_, u_, d_, curr_step));
    
    // The BS pricer
    auto last_step = [&](double S_T, double t)->double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);
        
        return option.Put();
    };
    
    // Get BS prices
    std::vector<double> V(this->GeneratePayoff(S[0], last_step, 0.));   // Again, the last param is a dummy.
    
    // Get early exercise payoffs
    std::vector<double> exercise_payoff(this->GeneratePayoff(S[0], payoff, curr_time));

    // Compare
    this->EarlyExUpdate(V, exercise_payoff);
    
    S[0].pop_back();
    S[0].pop_front();
    curr_step--;
    curr_time -= dt_;
    
    while (curr_step >= 2) {
        this->BacktrackEE(V, S, curr_step, payoff, curr_time);
        curr_step--;
        curr_time -= dt_;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];
    assert(V.size() == 3);
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_step--;
    curr_time -= dt_;
    
    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_step--;
    curr_time -= dt_;
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V20, V21, V22);
}

TreeResult BTree::PIBSR(const std::function<double (double, double)>& payoff) const {
    TreeResult res1 = this->PIBS(payoff);
    
    BTree shorter_tree(S0_, sigma_, T_, steps_ / 2, r_, q_);
    
    TreeResult res2 = shorter_tree.PIBS(payoff);
    
    double value = res1.value * 2 - res2.value;
    double delta = res1.delta * 2 - res2.delta;
    double gamma = res1.gamma * 2 - res2.gamma;
    double theta = res1.theta * 2 - res2.theta;
    
    return TreeResult({value, delta, gamma, theta});
}

TreeResult BTree::EEBSR(const std::function<double (double, double)>& payoff) const {
    TreeResult res1 = this->EEBS(payoff);
    
    BTree shorter_tree(S0_, sigma_, T_, steps_ / 2, r_, q_);
    
    TreeResult res2 = shorter_tree.EEBS(payoff);
    
    double value = res1.value * 2 - res2.value;
    double delta = res1.delta * 2 - res2.delta;
    double gamma = res1.gamma * 2 - res2.gamma;
    double theta = res1.theta * 2 - res2.theta;
    
    return TreeResult({value, delta, gamma, theta});
}

TreeResult BTree::BarrierVanilla(const std::function<double (double, double)>& payoff, const Barrier& barrier, double B) const {
    // TODO: Only down-and-out
    assert(barrier == DownAndOut);
    
    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::array<std::deque<double>, 2> S(this->GenerateSMeshEE(S0_, u_, d_, steps_));
    
    double curr_time = T_;
    std::size_t curr_step = steps_;
    std::vector<double> V(this->GeneratePayoffBarrier(S[0], payoff, curr_time, barrier, B));
    S[0].pop_back();
    S[0].pop_front();
    curr_step--;
    curr_time -= dt_;
    
    while (curr_step >= 2) {
        // During backtracking, curr_step and curr_time are exactly the current step (0-indexed) and the current time
        this->BacktrackBarrier(V, S, curr_step, barrier, B);
        curr_step--;
        curr_time -= dt_;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];
    assert(V.size() == 3);
    
    // Backtrack
    this->BacktrackBarrier(V, S, curr_step, barrier, B);
    curr_step--;
    curr_time -= dt_;
    
    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    
    // Backtrack
    this->BacktrackBarrier(V, S, curr_step, barrier, B);
    curr_step--;
    curr_time -= dt_;
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V20, V21, V22);
}

#endif /* BTree_cpp */

#include "TTree.hpp"
#include <cmath>
#include<cassert>
#include <algorithm>
#include <iostream>

TTree::TTree(double S0, double sigma, double T, std::size_t steps, double r, double q) : S0_(S0), sigma_(sigma), T_(T), steps_(steps), r_(r), q_(q),
    dt_(T / steps),
    u_(std::exp(sigma * std::sqrt(3. * dt_))),
    d_(1. / u_),
    r_disc_(std::exp(-r * dt_))
{
//    double erq = std::exp((r_ - q_) * dt_ / 2.);
//    double esig = std::exp(sigma_ * std::sqrt(dt_ / 2.));
//    pu_ = (erq - 1. / esig) / (esig - 1. / esig);
//    pu_ *= pu_;
//    pd_ = (esig - erq) / (esig - 1. / esig);
//    pd_ *= pd_;
    
    double temp = (r_ - q_ - sigma_ * sigma_ * .5) * std::sqrt(dt_ / (12. * sigma_ * sigma_));
    pu_ = 1. / 6. + temp;
    pd_ = 1. / 6. - temp;
    
    disc_pu = pu_ * r_disc_;
    disc_pd = pd_ * r_disc_;
    disc_pm = 2. / 3. * r_disc_;
//    disc_pm = (1. - pu_ - pd_) * r_disc_;
}

TTree::TTree(const EuropeanOption& option, std::size_t steps) : TTree(option.S_, option.sigma_, option.T_, steps, option.r_, option.q_) {}

TreeResult TTree::PathIndependentOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const {
    switch (modifier) {
        case TreeModifier::vanilla:
            // Vanilla tree
            return this->PIVanilla(payoff);
            break;
//        case TreeModifier::average:
//            // Average tree
//            return this->PIAvg(payoff);
//            break;
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

TreeResult TTree::EarlyExerciseOption(const std::function<double (double, double)>& payoff, const TreeModifier& modifier) const {
    switch (modifier) {
        case TreeModifier::vanilla:
            // Vanilla tree
            return this->EEVanilla(payoff);
            break;
//        case TreeModifier::average:
//            // Average tree
//            return this->EEAvg(payoff);
//            break;
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

std::deque<double> TTree::GenerateSMesh(double S0, double u, double d, std::size_t steps) const {
    
    // Generate asset price at maturity
    std::deque<double> S({S0 * std::pow(u, static_cast<double>(steps))});
    
    for (std::size_t i = 0; i < 2 * steps; i++) {
        S.push_back(S.back() * d);
    }
    
    return S;
}

std::vector<double> TTree::GeneratePayoff(const std::deque<double>& S, const std::function<double (double, double)>& payoff, double t) const {
    auto payoff_t = std::bind(payoff, std::placeholders::_1, t);
    
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), payoff_t);
    
    return V;
}

std::vector<double> TTree::GeneratePayoff(const std::deque<double> &S, const std::function<double (double, double)> &payoff, double t, const Barrier& barrier, double B) const {
    
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

void TTree::EarlyExUpdate(std::vector<double>& V, const std::vector<double>& earlyex_payoff) const {
    assert(V.size() == earlyex_payoff.size());
    
    auto earlyex_it = earlyex_payoff.cbegin();
    for (auto Vit = V.begin(); Vit != V.end(); Vit++) {
        *Vit = std::max(*Vit, *earlyex_it);
        earlyex_it++;
    }
}

void TTree::BacktrackPI(std::vector<double>& V) const {
    std::vector<double> new_V;
    
    for (auto Vit = V.cbegin() + 1; Vit + 1 != V.cend(); Vit++) {
        // Get value of the last node by taking the expectation and discounting
        // V = disc(p_u * V_u + p_m * V_m + p_d * V_d)
        new_V.push_back(disc_pu * (*(Vit - 1)) + disc_pm * (*Vit) + disc_pd * (*(Vit + 1)));
    }
    
    // Replace V
    V = std::move(new_V);
}

void TTree::BacktrackEE(std::vector<double>& V, std::deque<double>& S, size_t curr_step, const std::function<double (double, double)>& payoff, double curr_time) const {
    // Backtrack as if European
    this->BacktrackPI(V);
    
    // Find the payoff if we exercise the option at this step
    std::vector<double> exercise_payoff(this->GeneratePayoff(S, payoff, curr_time));
    
    this->EarlyExUpdate(V, exercise_payoff);
    
    if (S.size() >= 2) {
        // Shrink the current S since we will have fewer outcomes in the next iteration.
        S.pop_front();
        S.pop_back();
    }
    
}

void TTree::BacktrackBarrier(std::vector<double>& V, std::deque<double>& S, const Barrier& barrier, double B) const {
    
    assert(barrier == DownAndOut);
    // Backtrack as if European
    this->BacktrackPI(V);
    
    // Set to zero if S is down
    std::transform(S.cbegin(), S.cend(), V.cbegin(), V.begin(), [&](double s, double v)->double {
        if (s <= B) {
            return 0.;
        } else {
            return v;
        }
    });
    
    if (S.size() >= 2) {
        // Shrink the current S since we will have fewer outcomes in the next iteration.
        S.pop_front();
        S.pop_back();
    }
    
}

TreeResult TTree::GenTreeResult(double V00, double V10, double V11, double V12, double V20, double V22, double V24) const {
    
    double delta = (V10 - V12) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V22) / (S0_ * (u_ * u_ - 1.)) - ((V22 - V24) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V11 - V00)/ dt_;
    
    return TreeResult({V00, delta, gamma, theta});
}

TreeResult TTree::PIVanilla(const std::function<double (double, double)>& payoff) const {
    
    // Generate asset price at maturity
    std::deque<double> S(this->GenerateSMesh(S0_, u_, d_, steps_));
    
    // Generate derivative value at maturity
    std::vector<double> V(this->GeneratePayoff(S, payoff, T_));
    
    // Backtrack until the end
    while (V.size() > 5) {
        this->BacktrackPI(V);
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
//    double V21 = V[1];
    double V22 = V[2];
//    double V23 = V[3];
    double V24 = V[4];
    assert(V.size() == 5);

    // Backtrack
    this->BacktrackPI(V);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    this->BacktrackPI(V);
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V12, V20, V22, V24);
}

TreeResult TTree::EEVanilla(const std::function<double (double, double)>& payoff) const {
    // Generate asset price at maturity
    std::deque<double> S(this->GenerateSMesh(S0_, u_, d_, steps_));
    
    // Generate derivative value at maturity
    double curr_time = T_;
    std::size_t curr_step = steps_;
    std::vector<double> V(this->GeneratePayoff(S, payoff, curr_time));
    
    S.pop_front();
    S.pop_back();
    
    curr_time -= dt_;
    curr_step--;
    
    // Backtrack until the end
    while (curr_step >= 2) {
        this->BacktrackEE(V, S, curr_step, payoff, curr_time);
        curr_time -= dt_;
        curr_step--;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
//    double V21 = V[1];
    double V22 = V[2];
//    double V23 = V[3];
    double V24 = V[4];
    assert(V.size() == 5);
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_time -= dt_;
    curr_step--;
    
    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_time -= dt_;
    curr_step--;
    
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V12, V20, V22, V24);
}

TreeResult TTree::PIBS(const std::function<double (double, double)>& payoff) const {
    
    std::size_t curr_step = steps_;
    curr_step--;    // Don't need the first step
    
    // Generate asset price at maturity
    std::deque<double> S(this->GenerateSMesh(S0_, u_, d_, curr_step));
    
    // Generate the BS pricer
    auto last_step = [&](double S_T, double t)->double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);
        
        return option.Put();
    };
    
    // Generate derivative value at (maturity - dt)
    std::vector<double> V(this->GeneratePayoff(S, last_step, 0.));  // Last param is dummy
    
    curr_step--;
    
    // Backtrack until the end
    while (curr_step >= 2) {
        this->BacktrackPI(V);
        curr_step--;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
//    double V21 = V[1];
    double V22 = V[2];
//    double V23 = V[3];
    double V24 = V[4];
    assert(V.size() == 5);
    
    // Backtrack
    this->BacktrackPI(V);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    this->BacktrackPI(V);
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V12, V20, V22, V24);
}

TreeResult TTree::EEBS(const std::function<double (double, double)>& payoff) const {
    std::size_t curr_step = steps_;
    double curr_time = T_;
    curr_step--;    // We do not need the last step here
    curr_time -= dt_;
    
    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::deque<double> S(this->GenerateSMesh(S0_, u_, d_, curr_step));
    
    // The BS pricer
    auto last_step = [&](double S_T, double t)->double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);
        
        return option.Put();
    };
    
    // Get BS prices
    std::vector<double> V(this->GeneratePayoff(S, last_step, 0.));   // Again, the last param is a dummy.
    
    // Get early exercise payoffs
    std::vector<double> exercise_payoff(this->GeneratePayoff(S, payoff, curr_time));

    // Compare
    this->EarlyExUpdate(V, exercise_payoff);
    
    S.pop_back();
    S.pop_front();
    curr_step--;
    curr_time -= dt_;
    
    while (curr_step >= 2) {
        this->BacktrackEE(V, S, curr_step, payoff, curr_time);
        curr_step--;
        curr_time -= dt_;
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
//    double V21 = V[1];
    double V22 = V[2];
//    double V23 = V[3];
    double V24 = V[4];
    assert(V.size() == 5);
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_step--;
    curr_time -= dt_;
    
    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];
    
    // Backtrack
    this->BacktrackEE(V, S, curr_step, payoff, curr_time);
    curr_step--;
    curr_time -= dt_;
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V12, V20, V22, V24);
}

TreeResult TTree::PIBSR(const std::function<double (double, double)>& payoff) const {
    TreeResult res1 = this->PIBS(payoff);
    
    TTree shorter_tree(S0_, sigma_, T_, steps_ / 2, r_, q_);
    
    TreeResult res2 = shorter_tree.PIBS(payoff);
    
    double value = res1.value * 2 - res2.value;
    double delta = res1.delta * 2 - res2.delta;
    double gamma = res1.gamma * 2 - res2.gamma;
    double theta = res1.theta * 2 - res2.theta;
    
    return TreeResult({value, delta, gamma, theta});
}

TreeResult TTree::EEBSR(const std::function<double (double, double)>& payoff) const {
    TreeResult res1 = this->EEBS(payoff);
    
    TTree shorter_tree(S0_, sigma_, T_, steps_ / 2, r_, q_);
    
    TreeResult res2 = shorter_tree.EEBS(payoff);
    
    double value = res1.value * 2 - res2.value;
    double delta = res1.delta * 2 - res2.delta;
    double gamma = res1.gamma * 2 - res2.gamma;
    double theta = res1.theta * 2 - res2.theta;
    
    return TreeResult({value, delta, gamma, theta});
}

TreeResult TTree::BarrierVanilla(const std::function<double (double, double)>& payoff, const Barrier& barrier, double B) const {
    
    // TODO: Only down-and-out
    assert(barrier == DownAndOut);
    
    // Generate asset price at maturity
    std::deque<double> S(this->GenerateSMesh(S0_, u_, d_, steps_));
    
    // Generate derivative value at maturity
    std::vector<double> V(this->GeneratePayoff(S, payoff, T_, barrier, B));
    
    S.pop_front();
    S.pop_back();
    
    // Backtrack until the end
    while (V.size() > 5) {
        this->BacktrackBarrier(V, S, barrier, B);
    }
    
    // Store nodes of Step 2
    double V20 = V[0];
//    double V21 = V[1];
    double V22 = V[2];
//    double V23 = V[3];
    double V24 = V[4];
    assert(V.size() == 5);

    // Backtrack
    this->BacktrackBarrier(V, S, barrier, B);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    this->BacktrackBarrier(V, S, barrier, B);
    double V00 = V[0];
    
    return this->GenTreeResult(V00, V10, V11, V12, V20, V22, V24);

    
}


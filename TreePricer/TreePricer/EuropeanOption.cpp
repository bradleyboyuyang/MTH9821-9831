//
//  EuropeanOption.cpp
//  BinomialTree
//
//  Created by Mingsen Wang on 2022/6/7.
//

#include "EuropeanOption.hpp"
#include <cmath>
#include <numbers>
//#include <iostream>
#include <algorithm>

double EuropeanOption::Phi(double t) const {
    return std::erfc(-t / std::sqrt(2.)) / 2.;
}

double EuropeanOption::Z(double t) const {
    return std::exp(-t * t / 2.) / std::sqrt(2. * std::numbers::pi);
}

EuropeanOption::EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q) : t_(t), S_(S), K_(K), T_(T), sigma_(sigma), r_(r), q_(q) {
    d1_ = (log(S/K) + (r - q + sigma * sigma / 2.) * (T - t)) / (sigma * std::sqrt(T-t));
    d2_ = d1_ - sigma * std::sqrt(T-t);
    
    Zd1_ = this->Z(d1_);
    Zd2_ = this->Z(d2_);
    Nd1_ = this->Phi(d1_);
    Nd2_ = this->Phi(d2_);
    
    q_disc_ = std::exp(-q * (T-t));
    r_disc_ = std::exp(-r * (T-t));
}

double EuropeanOption::Call() const {
    return S_ * q_disc_ * Nd1_ - K_ * r_disc_ * Nd2_;
}

double EuropeanOption::Put() const {
    return -S_ * q_disc_ * (1. - Nd1_) + K_ * r_disc_ * (1. - Nd2_);
}

double EuropeanOption::DeltaCall() const {
    return q_disc_ * Nd1_;
}
double EuropeanOption::DeltaPut() const {
    return -q_disc_ * (1. - Nd1_);
}

double EuropeanOption::GammaCall() const {
    return q_disc_ / (S_ * sigma_ * std::sqrt(T_ - t_)) * Zd1_;
}
double EuropeanOption::GammaPut() const {
    return this->GammaCall();
}

double EuropeanOption::ThetaCall() const {
    double res = -(S_ * sigma_ * q_disc_) / (2. * std::sqrt(T_ - t_)) * Zd1_;
    
    res += q_ * S_ * q_disc_ * Nd1_;
    
    res += -r_ * K_ * r_disc_ * Nd2_;
    
    return res;
}

double EuropeanOption::ThetaPut() const {
    double res = -(S_ * sigma_ * q_disc_) / (2. * std::sqrt(T_ - t_)) * Zd1_;
    
    res += -q_ * S_ * q_disc_ * (1 - Nd1_);
    
    res += r_ * K_ * r_disc_ * (1 - Nd2_);
    
    return res;
}

double EuropeanOption::VegaCall() const {
    return S_ * q_disc_ * std::sqrt(T_ - t_) * Zd1_;
}

double EuropeanOption::VegaPut() const {
    return this->VegaCall();
}

double EuropeanOption::RhoCall() const {
    return K_ * (T_ - t_) * r_disc_ * Nd2_;
}

double EuropeanOption::RhoPut() const {
    return -K_ * (T_ - t_) * r_disc_ * (1. - Nd2_);
}

std::function<double(double, double)> EuropeanOption::CallPayoff() const {
    return [&](double S, double t)->double {
        return std::max(S - K_, 0.);
    };
}

std::function<double(double, double)> EuropeanOption::PutPayoff() const {
    return [&](double S, double t)->double {
        return std::max(K_ - S, 0.);
    };
}

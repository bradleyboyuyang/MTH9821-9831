//
//  main.cpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/11/22.
//

#include <iostream>
#include <iomanip>

#include "EuropeanOption.hpp"
#include "finite_difference_pricer.hpp"
#include "finite_difference_pricer_discrete_dividend.hpp"

void PrintVector(const std::vector<double>& vec) {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}

void Euro_explicit();
void American_explicit();
//void American_boundary();

void Euro_implicit(double alpha = .45) {
    EuropeanOption option(0., 37., 40, .75, .28, .04, .02);
    FiniteDifferencePricer::alpha_temp_ = alpha;
    FiniteDifferencePricer pricer(37., 40, .75, .28, .04, .02);
    
    std::cout << std::setprecision(8);
    std::cout << "Backward Euler, alpha = " << alpha << std::endl;
    std::cout << option.Put() << std::endl;
    
    for (std::size_t i = 4; i < 5; i <<= 2) {
        PrintVector(pricer.Price(i, Implicit, European));
    }
}

void Euro_imex(double alpha = .45) {
    EuropeanOption option(0., 37., 40, .75, .28, .04, .02);
    FiniteDifferencePricer::alpha_temp_ = alpha;
    FiniteDifferencePricer pricer(37., 40, .75, .28, .04, .02);
    
    std::cout << std::setprecision(8);
    std::cout << "Crank-Nicolson, alpha = " << alpha << std::endl;
//    std::cout << option.Put() << std::endl;
    
    for (std::size_t i = 4; i < 5; i <<= 2) {
        PrintVector(pricer.Price(i, ImEx, European));
    }
}

void American_imex(double alpha = .45) {
    EuropeanOption option(0., 37., 40, .75, .28, .04, .02);
    FiniteDifferencePricer::alpha_temp_ = alpha;
    FiniteDifferencePricer pricer(37., 40, .75, .28, .04, .02);
    
    std::cout << std::setprecision(10);
    std::cout << "Crank-Nicolson, alpha = " << alpha << std::endl;
//    std::cout << option.Put() << std::endl;
    
    for (std::size_t i = 4; i < 5; i <<= 2) {
        PrintVector(pricer.Price(i, ImEx, American));
    }
}


void American_imex2(double alpha = .45) {
    EuropeanOption option(0., 42., 40., 1., .32, .04, .02);
    FiniteDifferencePricer::alpha_temp_ = alpha;
    FiniteDifferencePricer pricer(42., 40., 1., .32, .04, .02);
    
    std::cout << std::setprecision(10);
    std::cout << "Crank-Nicolson, alpha = " << alpha << std::endl;
//    std::cout << option.Put() << std::endl;
    
    for (std::size_t i = 4; i < 5; i <<= 2) {
        PrintVector(pricer.Price(i, ImEx, American));
    }
}

void AmericanDiscDiv(double alpha = .4) {
    FiniteDifferencePricerDiscreteDividend::alpha_temp_ = 4;
    
    EuropeanOption option(0., 52., 50., 1., .2, .03, 0.);
    FiniteDifferencePricerDiscreteDividend pricer(option, 5. / 12., .01);
    
    PrintVector(pricer.PriceCall(4, ImEx, American));
}

void Final() {
    EuropeanOption option(0., 42., 45., 8. / 12., 0.18820914, .04, .02);
    FiniteDifferencePricer pricer(option);
    
    std::cout << pricer.Price(200, 40, ImEx, American)[0] << std::endl;
}

double IV() {
//    double S = 42.;
//    double K = 45.;
//    double T = 8. / 12.;
//    double r = 0.04;
//    double q = 0.02;
    
    double V0 = 4.09;
    
    double tol_approx = std::pow(10., -4);
    
    double sig_last = .1;
    double sig_new = .5;
    
    EuropeanOption option_new(0., 42., 45., 8. / 12., sig_new, .04, .02);
    FiniteDifferencePricer pricer_new(option_new);
    double v_new = pricer_new.Price(2000, 400, ImEx, American).front();
    
    EuropeanOption option_last(0., 42., 45., 8. / 12., sig_last, .04, .02);
    FiniteDifferencePricer pricer_last(option_last);
    double v_last = pricer_last.Price(2000, 400, ImEx, American).front();
    
    while (std::abs(sig_new - sig_last) >= tol_approx) {
        // Secant method
        double sig_newest = sig_new - (v_new - V0) * (sig_new - sig_last) / (v_new - v_last);

        sig_last = sig_new;
        sig_new = sig_newest;
        v_last = v_new;
        EuropeanOption option_newest(0., 42., 45., 8. / 12., sig_newest, .04, .02);
        FiniteDifferencePricer pricer_newest(option_newest);
        v_new = pricer_newest.Price(2000, 400, ImEx, American).front();
        std::cout << "Sig:\t" << sig_new << std::endl;
    }
    
    return sig_new;
//    EuropeanPut new_put(S_, K_, T_, sig_new, r_, q_);
//    EuropeanPut last_put(S_, K_, T_, sig_last, r_, q_);
//
//    double tol_approx = std::pow(10, -4);
//
//    while (std::abs(sig_new - sig_last) >= tol_approx) {
//        EuropeanOption option_new(0., 42., 45., 8. / 12., .3, .04, .02);
//        FiniteDifferencePricer pricer_new(option);
//
//        EuropeanPut new_put(S_, K_, T_, sig_new, r_, q_);
//        EuropeanPut last_put(S_, K_, T_, sig_last, r_, q_);
//
//        double v_new = new_put.BinomialTree(steps, modifier).value;
//        double v_last = last_put.BinomialTree(steps, modifier).value;
//
//        // Secant method
//        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);
//
//        sig_last = sig_new;
//        sig_new = sig_newest;
//    }
//
//    return sig_new;
}


int main(int argc, const char * argv[]) {
    
//    Euro_explicit();
//    American_explicit();
//    American_boundary();
//    Euro_implicit(.45);
//    Euro_imex(.45);
//    Euro_implicit(5);
//    Euro_imex(5);
    
//    American_imex(.45);
//    American_imex(5.);
    std::cout << std::fixed << std::setprecision(8);
    
    Final();
//    std::cout << IV() << std::endl;
    
    
    return 0;
}

// HW7
void Euro_explicit() {
    EuropeanOption option(0., 37., 40, .75, .28, .03, .015);
    FiniteDifferencePricer pricer(37., 40, .75, .28, .03, .015);
    
    std::cout << std::setprecision(8);
    std::cout << option.Put() << std::endl;
    for (std::size_t i = 4; i < 257; i <<= 2) {
        PrintVector(pricer.Price(i, Explicit, European));
    }
}

void American_explicit() {
    EuropeanOption option(0., 42., 40, .75, .32, .02, .04);
    FiniteDifferencePricer pricer(option);
    
    std::cout << std::setprecision(8);
    
    for (std::size_t i = 16; i < 17; i <<= 2) {
        PrintVector(pricer.Price(i, Explicit, American));
    }
}
//
//void American_boundary() {
//    FiniteDifferencePricer pricer(42., 40, .75, .32, .04, .02);
//
//    std::cout << std::setprecision(4);
//
//    std::vector<double> t;
//    std::vector<double> Sopt;
//    std::tie(t, Sopt) = pricer.AmericanPut_EarlyExDomain(16);
//    PrintVector(t);
//    PrintVector(Sopt);
//}

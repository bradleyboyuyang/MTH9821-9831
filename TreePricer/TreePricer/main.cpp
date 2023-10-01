//
//  main.cpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#include <iostream>
#include "EuropeanPut.hpp"
#include "AmericanPut.hpp"
#include <iomanip>
#include <cmath>
#include "HW11.hpp"

//void TestBasicBTree() {
//    EuropeanOption option(0., 54., 50., 1., .29, .0375, .01);
//    EuropeanPut europut(option);
//    AmericanPut ameriput(option);
//
//    std::cout << "European" << std::endl;
//    europut.BinomialTree(1280, vanilla).Print();
//    europut.BinomialTree(1280, average).Print();
//    europut.BinomialTree(1280, BS).Print();
//    europut.BinomialTree(1280, BSR).Print();
//
//    std::cout << "American" << std::endl;
//    ameriput.BinomialTree(1280, vanilla, false).Print();
//    ameriput.BinomialTree(1280, average, false).Print();
//    ameriput.BinomialTree(1280, BS, false).Print();
//    ameriput.BinomialTree(1280, BSR, false).Print();
//
//    std::cout << "American, variance reduction" << std::endl;
//    ameriput.BinomialTree(1280, vanilla, true).Print();
//    ameriput.BinomialTree(1280, average, true).Print();
//    ameriput.BinomialTree(1280, BS, true).Print();
//    ameriput.BinomialTree(1280, BSR, true).Print();
//}
//
//void TestBasicTTree() {
//    EuropeanOption option(0., 41., 39., 1., .25, .03, .005);
//    EuropeanPut europut(option);
//    AmericanPut ameriput(option);
//
//    std::cout << "European" << std::endl;
//    europut.TrinomialTree(1280, vanilla).Print();
//    europut.TrinomialTree(1280, BS).Print();
//    europut.TrinomialTree(1280, BSR).Print();
//
//    std::cout << "American" << std::endl;
//    ameriput.TrinomialTree(1280, vanilla, false).Print();
//    ameriput.TrinomialTree(1280, BS, false).Print();
//    ameriput.TrinomialTree(1280, BSR, false).Print();
//
//    std::cout << "American, variance reduction" << std::endl;
//    ameriput.TrinomialTree(1280, vanilla, true).Print();
//    ameriput.TrinomialTree(1280, BS, true).Print();
//    ameriput.TrinomialTree(1280, BSR, true).Print();
//}
//
//void TestDiscreteDividendEuropeanPut() {
//    // 2% after 2 months
//    EuropeanOption option(0., 50., 55.55, .25, .2, .02, 0.);
//    EuropeanPut europut(option);
//    Dividend proportional({{2./12.}, {.02}});
//    Dividend fixed;
//    europut.BinomialTree(1536, vanilla, proportional, fixed).Print();
//
//    EuropeanOption option2(0., 49., 55.55, .25, .2, .02, 0.);
//    EuropeanPut europut2(option2);
//    std::cout << option2.Put() << '\t' << option2.DeltaPut() << '\t' << option2.GammaPut() << '\t' << option2.ThetaPut() << std::endl;
//    europut2.BinomialTree(1536, vanilla).Print();
//}
//
//void TestDiscreteDividendAmericanPut() {
//    // 2% after 2 months
//    EuropeanOption option(0., 50., 55.55, .25, .2, .02, 0.);
//    AmericanPut europut(option);
//    Dividend proportional({{2./12.}, {.02}});
//    Dividend fixed;
//    europut.BinomialTree(6, vanilla, proportional, fixed).Print();
//}

double DaO(const EuropeanOption& option, double B) {
    double a = (option.r_ - option.q_) / (option.sigma_ * option.sigma_) - .5;
    
    EuropeanOption option2(0., B * B / option.S_, option.K_, option.T_, option.sigma_, option.r_, option.q_);
    
    double V = option.Call() - std::pow((B / option.S_), 2. * a) * option2.Call();
    
    return V;
}

int FindOptimalN(const EuropeanOption& option, int k, double B) {
    return std::floor(k * k * 3. * option.sigma_ * option.sigma_ * option.T_ / std::pow(std::log(option.S_ / B), 2));
}

void HW12() {
    EuropeanOption option(0., 50., 48., 10. / 12., .30, .02, .01);
    EuropeanPut put(option);
    
    double B = 35.;
    
    double theoretical_val = DaO(option, B);
    std::cout << "Theoretical value: " << theoretical_val << std::endl;
    
    for (std::size_t N = 10; N <= 1000; N++) {
        double val = put.TrinomialTree(N, vanilla, DownAndOut, B).value;
        std::cout << N << '\t' << std::abs(val - theoretical_val) << std::endl;
    }
    
    for (int k = 2; k < 12; k++) {
        std::size_t optimal_N = FindOptimalN(option, k, B);
        double val = put.TrinomialTree(optimal_N, vanilla, DownAndOut, B).value;
        std::cout << optimal_N << '\t' << std::abs(val - theoretical_val) << std::endl;
    }
}

void Final() {
    EuropeanOption option(0., 50., 48., 10. / 12., .30, .02, .01);
    EuropeanPut put(option);
    
    double B = 45.;
    
    double theoretical_val = DaO(option, B);
    std::cout << "Theoretical value: " << theoretical_val << std::endl;
    
    for (std::size_t N = 10; N <= 15; N++) {
        double val = put.BinomialTree(N, vanilla, DownAndOut, B).value;
        double val2 = put.TrinomialTree(N, vanilla, DownAndOut, B).value;
        std::cout << val << '\t' << std::abs(val - theoretical_val) << '\t' << val2 << '\t' << std::abs(val2 - theoretical_val) << std::endl;
    }
}

int main(int argc, const char * argv[]) {
    
    std::cout << std::fixed << std::setprecision(8);
//    TestBasicBTree();
//    TestBasicTTree();
    
//    TestDiscreteDividendEuropeanPut();
//    TestDiscreteDividendAmericanPut();
    
//    HW11();
//    HW12();
    Final();
    
    return 0;
}

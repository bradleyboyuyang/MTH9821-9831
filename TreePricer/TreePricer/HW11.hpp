//
//  HW11.hpp
//  TreePricer
//
//  Created by 王明森 on 11/23/22.
//

#ifndef HW11_hpp
#define HW11_hpp

#include "EuropeanPut.hpp"
#include "AmericanPut.hpp"
#include <cmath>
#include <iostream>

void HW11_1a(double tol, std::size_t init_step) {
    EuropeanOption option(0., 50., 55.55, .25, .2, .02, 0.);
    EuropeanPut europut(option);
    AmericanPut ameriput(option);
    
    Dividend proportional({{2. / 12.}, {.02}});
    Dividend fixed;
    
    TreeResult res;
    std::size_t final_step;
    
    std::tie(res, final_step) = europut.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
    
    std::tie(res, final_step) = ameriput.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
}

void HW11_1b(double tol, std::size_t init_step) {
    EuropeanOption option(0., 50., 55.55, 7. / 12., .2, .02, 0.);
    EuropeanPut europut(option);
    AmericanPut ameriput(option);
    
    Dividend proportional({{2. / 12., 4. / 12., 6. / 12.}, {.02, .02, .02}});
    Dividend fixed;
    
    TreeResult res;
    std::size_t final_step;
    
    std::tie(res, final_step) = europut.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
    
    std::tie(res, final_step) = ameriput.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
}

void HW11_2a(double tol, std::size_t init_step) {
    EuropeanOption option(0., 50., 55.55, .25, .2, .02, 0.);
    EuropeanPut europut(option);
    AmericanPut ameriput(option);
    
    Dividend proportional;
    Dividend fixed({{2. / 12.}, {.75}});
    
    TreeResult res;
    std::size_t final_step;
    
    std::tie(res, final_step) = europut.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
    
    std::tie(res, final_step) = ameriput.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
}

void HW11_2b(double tol, std::size_t init_step) {
    EuropeanOption option(0., 50., 55.55, 7. / 12., .2, .02, 0.);
    EuropeanPut europut(option);
    AmericanPut ameriput(option);
    
    Dividend proportional;
    Dividend fixed({{2. / 12., 4. / 12., 6. / 12.}, {.75, .75, .75}});
    
    TreeResult res;
    std::size_t final_step;
    
    std::tie(res, final_step) = europut.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
    
    std::tie(res, final_step) = ameriput.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
}

void HW11_3(double tol, std::size_t init_step) {
    EuropeanOption option(0., 50., 55.55, 7. / 12., .2, .02, 0.);
    EuropeanPut europut(option);
    AmericanPut ameriput(option);
    
    Dividend proportional({{4. / 12.}, {.02}});
    Dividend fixed({{2. / 12., 6. / 12.}, {.75, .25}});
    
    TreeResult res;
    std::size_t final_step;
    
    std::tie(res, final_step) = europut.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
    
    std::tie(res, final_step) = ameriput.BinomialTree(init_step, vanilla, proportional, fixed, tol);
    std::cout << res.value << '\t' << final_step << '\t' << res.delta << '\t' << res.gamma << '\t'<< res.theta << std::endl;
}



void HW11() {
    double tol = std::pow(10., -4);
    HW11_1a(tol, 6);
    HW11_1b(tol, 7);
    HW11_2a(tol, 6);
    HW11_2b(tol, 7);
    HW11_3(tol, 7);
    
    
}

#endif /* HW11_hpp */

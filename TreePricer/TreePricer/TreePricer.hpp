//
//  TreePricer.hpp
//  TreePricer
//
//  Created by 王明森 on 11/22/22.
//

#ifndef TreePricer_hpp
#define TreePricer_hpp

#include <iostream>
#include <vector>

struct TreeResult {
    double value;
    double delta;
    double gamma;
    double theta;
    
    void Print() const {
        std::cout << value << '\t' << delta << '\t' << gamma << '\t' << theta << std::endl;
    }
};

enum TreeModifier {
    vanilla,
    average,    // Average tree
    BS,         // Black-Scholes
    BSR         // Black-Scholes with Richardson extrapolation
};

struct Dividend {
    std::vector<double> dates;
    std::vector<double> dividends;
    
    void Print() const {
        std::cout << "t: ";
        for (double t : dates) {
            std::cout << t << ' ';
        }
        std::cout << "\nd: ";
        for (double d : dividends) {
            std::cout << d << ' ';
        }
        std::cout << std::endl;
    }
};

//void PrintVector(const std::vector<double>& vec) {
//    for (double elem : vec) {
//        std::cout << elem << '\t';
//    }
//    std::cout << std::endl;
//}

enum Barrier {
    UpAndIn,
    UpAndOut,
    DownAndIn,
    DownAndOut,
};

#endif /* TreePricer_hpp */

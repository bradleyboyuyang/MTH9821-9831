#include<iostream>
#include"bsm.hpp"

int main(){
    double sum =0;
    double var = 0;;
    for(int i = 1; i < 100000; i++){
        sum +=bsm(i/100000.);
        
        std::cout << i << "   " << bsm(i/100000.) << std::endl;
    }
    std::cout << sum/100000. << std::endl;
    for(int i = 1; i < 100000; i++){
        var +=( bsm(i/100000.) - sum/100000.)*( bsm(i/100000.) - sum/100000.);
    }
    std::cout << var/100000. << std::endl;
    return 0;
}
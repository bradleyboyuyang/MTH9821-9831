#ifndef BSM
#define BSM
#include<cmath>

double bsm(double u){
    double a_0 = 2.50662823884;     double b_0 = -8.47351093090;
    double a_1 = -18.61500062529;   double b_1 = 23.08336743743;
    double a_2 = 41.39119773534;    double b_2 = -21.06224101826;
    double a_3 = -25.44106049637;   double b_3 = 3.13082909833;

    double c_0 = 0.3374754822726147; double c_5 = 0.0003951896511919;
    double c_1 = 0.9761690190917186; double c_6 = 0.0000321767881768;
    double c_2 = 0.1607979714918209; double c_7 = 0.0000002888167364;
    double c_3 = 0.0276438810333863; double c_8 = 0.0000003960315187;
    double c_4 = 0.0038405729373609;
    double x;
    double r;
    double y = u -0.5;
    if(std::abs(y) < 0.42){
        r = y*y;
        x = y*(((a_3 * r + a_2)*r + a_1)*r + a_0)/ 
            ((((b_3 * r + b_2)*r + b_1)*r + b_0)*r + 1.0);
    }else{
        r = u;
        if(y > 0){
            r = 1 - u;
        }
        r = std::log(-std::log(r));
        x = c_0 + r*(c_1 + r*(c_2 + r *(c_3 + r * 
            c_4 + r * (c_5 + r*(c_6 + r * (c_7 + r*c_8))))));
        if(y < 0){
            x = -x;
        }
    }
        
    
    return x;

}
#endif
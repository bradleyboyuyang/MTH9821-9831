#include"FiniteDifferencesPricer.h"
#include <fstream>
#include <iomanip>

using namespace std;

int main() {
    DAO_Call_Euro_Option op(European, call, 42, 40, 7.0/12, 0.28, 0.015, 0.04,36);
    
    std::ofstream out;
    out.open("data_1.csv");
    
    //question 1
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 0.4, M_, CN_SOR);
        out << setprecision(10) << p.alpha << "," << p.x_left <<"," << p.x_right <<"," << p.N <<"," << p.delta_x <<"," << p.delta_tau <<endl;
    }
    out << endl << endl;
    
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 4, M_, CN_SOR);
        out << setprecision(10) << p.alpha << "," << p.x_left <<"," << p.x_right <<"," << p.N <<"," << p.delta_x <<"," << p.delta_tau <<endl;
    }
    out << endl << endl;
    
    
    //question 2 Finite difference solvers
    //Forward Euler, for alpha_temp = 0.4
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 0.4, M_, Forward);
        out << setprecision(10) << p.u_value << "," << p.v_approx <<","<< p.error_pointwise <<"," << p.delta_approx <<"," << p.gamma_approx <<"," << p.theta_approx <<endl;
    }
    out << endl << endl;
    
    //Backward Euler, for alpha_temp = 0.4
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 0.4, M_, Backward);
        out << setprecision(10) << p.u_value << "," << p.v_approx <<"," << p.error_pointwise <<"," << p.delta_approx <<"," << p.gamma_approx <<"," << p.theta_approx <<endl;
    }
    out << endl << endl;
    
    //Backward Euler, for alpha_temp = 4
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 4, M_, Backward);
        out << setprecision(10) << p.u_value << "," << p.v_approx <<"," << p.error_pointwise <<"," << p.delta_approx <<"," << p.gamma_approx <<"," << p.theta_approx <<endl;
    }
    out << endl << endl;
    
    
    //Crank-Nicolson, for alpha_temp = 0.4
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 0.4, M_, CN_SOR);
        out << setprecision(10) << p.u_value << "," << p.v_approx <<"," << p.error_pointwise <<"," << p.delta_approx <<"," << p.gamma_approx <<"," << p.theta_approx <<endl;
    }
    out << endl << endl;
    
    //Backward Euler, for alpha_temp = 4
    for (int M_ = 4;M_ <= 256;M_ *= 4) {
        cout << M_ << endl;
        FiniteDifferencePircer p(op, 4, M_, CN_SOR);
        out << setprecision(10) << p.u_value << "," << p.v_approx <<"," << p.error_pointwise <<"," << p.delta_approx <<"," << p.gamma_approx <<"," << p.theta_approx <<endl;
    }
    out << endl << endl;
    
    
    //question 3
    FiniteDifferencePircer p1(op, 0.4, 4, Forward);
    for (int m = 0;m <= p1.M;++m) {
        for (int n = 0;n <= p1.N;++n) {
            out << setprecision(10) << p1.Domin[n][m] << ",";
        }
        out<<endl;
    }
    
    out << endl << endl;
    
    FiniteDifferencePircer p2(op, 0.4, 4, Backward);
    for (int m = 0;m <= p2.M;++m) {
        for (int n = 0;n <= p2.N;++n) {
            out << setprecision(10) << p2.Domin[n][m] << ",";
        }
        out<<endl;
    }
    
    
    out.close();
    return 0;
}


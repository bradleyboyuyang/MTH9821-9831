//
// Created by Ziyuan Li on 2023/10/9.
//
#include "mc_regression.h"
#include "regression.h"
#include <limits>
#include "american_pricers.h"
#include <iostream>
namespace montecarlo {
    using namespace linalg::regression;
    double regression_pricer_backward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate,
                                      regression::MCRegression& mc_regression,
                                      const MonteCarloRegressionMethod& method){
        // method: Assume that this is set to Tsitsiklis VanRoy for now, we will discuss the Longstaff Schwartz case later.
        //w: The payoff of the option will be max(w(S − K), 0), so w = 1 means call, and w = −1 means put.
        //mc regression: An instance of PolynomialMCRegression or HermiteMCRegression. This object will estimate the continuation values Ci(S) and store the coefficients βij seen in equation (4).
        //paths: An array with N rows and M columns. Each of the N rows of paths represents a randomly generated path. More formally, its entries are:


        //In the body of regression pricer backward, we will need the following arrays, each of length N:
        //• E will be the immediate Exercise value
        //• C will be the estimated Continuation value, computed at each step by regression
        // • P will be the actual Payoff value
        // The Tsitsiklis–Van Roy algorithm can be summarised by three equations:
        // Et =max(w(St −K),0)
        // Ct = proj(P(t+δt) * e^(−rδt)|span(f1(St),f2(St),...fB(St))) (7)
        //  Pt= Et if t = T
        //  Pt= max(Et,Ct) if 0<t<T
        // The proj operator in equation 7 means “calculate the fitted values of a linear regression”. To compute this, please use the fit predict method of the mc_regression argument.
        // When t = 0, the stock price is known: S(0). In this case, the orthogonal projection makes C0 the average of Pδte−rδt across paths. In this case the computation should be done in the fit predict at 0 method.
        // The algorithm returns a single value, P0 := max(E0, C0).
        //Implement regression pricer backward with a for loop, using the recursive formulation above.
        // With a small code change, your function can also compute the Longstaff– Schwartz value. The only difference is equation (8) which should be modified to:
        //Pt = Et if t = T or if Et ≥Ct and0<t<T
        // Pt =  P(t+δt)*e^(−rδt) if Et <Ct and 0<t<T
        //Use the regression pricer backward function’s method argument to enable the caller to switch between the two algorithms.
        arr P;
        arr E;
        arr C;
        int M = int(paths.cols());
        double dt = maturity / M;
        if (method == MonteCarloRegressionMethod::Tsitsiklis_VanRoy){
            for (int j = M-1; j > 0; j--){
                E = w*(paths.col(j)-strike).max(0);
                if (j == M-1){
                    P = E;
                }
                else{
                    C = mc_regression.fit_predict(j, paths.col(j), P * exp(-interest_rate * dt));
                    P = E.max(C);}

            }
            double C0 = mc_regression.fit_predict_at_0(P * exp(-interest_rate * dt));
            double P0 = std::max(C0,w*(spot-strike));
            return P0;
        }

        else if (method == MonteCarloRegressionMethod::Longstaff_Schwartz){
            for (int j = M-1; j > 0; j--){
                E = w*(paths.col(j)-strike).max(0);
                if (j == M-1){
                     P = E;
                }
                else{
                    C = mc_regression.fit_predict(j, paths.col(j), P * exp(-interest_rate * dt));
                    // P = where(E>=C,E,P * exp(-interest_rate * dt));
                    for (int k = 0; k < E.size(); k++){
                        if (E(k) >= C(k)){
                            P(k) = E(k);
                        }
                        else{
                            P(k) = P(k) * exp(-interest_rate * dt);
                        }
                    }
                }
            }
            double C0 = mc_regression.fit_predict_at_0(P * exp(-interest_rate * dt));
            double P0 = std::max(C0,w*(spot-strike));
            return P0;
        }
        else{
            throw std::invalid_argument("Invalid method");
        }
    }


    double regression_pricer_forward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                 const double interest_rate,
                                 const regression::MCRegression& mc_regression){
        //    When running regression pricer backward, it not only prices the option, but also modifies the supplied MCRegression object2. After pricing, the regression object contains all the fitted coefficients, so all the information needed to calculate the continuation values Ci(S).
        //            This allows us to easily implement the “forward version” of our pricer:
        // As before we will use an array P of length N for the actual payoff values.
        // For each path k, this method starts at time t0 = 0, and step by step compares the immediate exercise value K − S to the estimated continuation value Ci(S).
        // If at time j*δt we find that K −S_{j*δt,k} exceeds Cj(S_{jδt,k}), the option is exercised, we set Pk = e^(−jδt)*r*(K − S), break the forward loop and move on to the next path, starting at time 0.
        //Return the average payoff value across paths: mean(P ).
        //Your implentation should use both MCRegression::predict and MCRegression:: predict at 0.

        int M = int(paths.cols());
        int N = int(paths.rows());
        double dt = maturity / M;
        arr P(N);
        arr2 C(N, M);
        double E;
        arr C0(N);
        C0.fill(mc_regression.predict_at_0());
        C.col(0) = C0;
        for(int j = 1; j < M-1 ; j++){
            C.col(j) = (mc_regression.predict(j, paths.col(j)));
        }

        C.col(M-1)  = w * (paths.col(M-1) - strike);

        for (int k = 0; k < N; k++){
            for (int j = 0; j < M; j++) {
                E = w * (paths(k, j) - strike);
                if (E > C(k, j)) {
                    P(k) = E * exp(-interest_rate * j * dt);
                    break;
                }
                P(k) = 0;
            }
        }
        return P.mean();
    }
}
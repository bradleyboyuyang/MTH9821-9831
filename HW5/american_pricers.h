#ifndef MC_AMERICAN_PRICERS_H
#define MC_AMERICAN_PRICERS_H

#include "mc_regression.h"
//#include "propre/types.h"
//#include "finance/instruments.h"

namespace montecarlo {
	typedef Eigen::ArrayXd arr;
	typedef Eigen::ArrayXXd arr2;

    enum MonteCarloRegressionMethod {Longstaff_Schwartz, Tsitsiklis_VanRoy};

    double regression_pricer_backward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate,
                                      regression::MCRegression& mc_regression,
                                      const MonteCarloRegressionMethod& method);
    double regression_pricer_forward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate,
                                      const regression::MCRegression& mc_regression);
}

#endif // MC_AMERICAN_PRICERS_H

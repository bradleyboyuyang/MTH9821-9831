#include "mc_regression.h"
#include "regression.h"
#include <limits>

namespace montecarlo {
	
    namespace regression {

        using namespace linalg::regression;

        MCRegression::MCRegression(const uint times):
            _times(times),
            _coeff(std::numeric_limits<double>::quiet_NaN()) {}

        MCRegression::~MCRegression() {}

        double MCRegression::fit_predict_at_0(const arr& y) {
            return _coeff = y.mean();
        }

        double MCRegression::predict_at_0() const {
            return _coeff;
        }

        PolynomialMCRegression::PolynomialMCRegression(const uint times, const uint degree, const bool precondition):
            MCRegression(times),
            degree(degree),
            precondition(precondition) {
            _coeffs.resize(times);
        }

        arr PolynomialMCRegression::fit_predict(const uint time, const arr& x, const arr& y) {
            if (time == 0) {
                throw std::out_of_range ("At time 0, fit_predict_at_0 should be used.");
            }
            arr result;
            std::tie(_coeffs.at(time), result) = fit_linear_regression(vandermonde(x, degree), y, precondition);
            return result;
        }


        arr PolynomialMCRegression::predict(const uint time, const arr& x) const {
            if (time == 0) {
                throw std::out_of_range ("At time 0, predict_at_0 should be used.");
            }
            return evaluate_polynomial(_coeffs.at(time), x);
        }

        HermiteMCRegression::HermiteMCRegression(const uint times, const uint degree, const bool precondition,
                                                       const bool standardize)
            /* Implement a constructor similar to that of PolynomialMCRegression */

        arr HermiteMCRegression::fit_predict(const uint time, const arr &x, const arr &y) {
			/* Implement a Hermite version of PolynomialMCRegression::fit_predict */
            if (time == 0) {
                throw std::out_of_range ("At time 0, fit_predict_at_0 should be used.");
            }
        }

        arr HermiteMCRegression::predict(const uint time, const arr& x) const {
            /* Implement a Hermite version of PolynomialMCRegression::predict */
            if (time == 0) {
                throw std::out_of_range ("At time 0, fit_predict_at_0 should be used.");
            }
        }
    }
}


#include "regression.h"
#include <tuple>

namespace linalg {

	namespace regression {

		using namespace Eigen;

		arr evaluate_polynomial(const arr& coefficients, const arr& x) {
			/* Evaluate a polynomial at an array of points
			 *
			 * coefficients(0) is the constant term, coefficients(1) is the linear term, etc.
			 * Use Horner's rule
			 * Code should be 5-10 lines
			 */
		}

		arr evaluate_hermite_polynomial(const arr& coefficients, const arr& x) {
			/* Evaluate a polynomial expressed in the Hermite basis at an array of points
			 *
			 * coefficients(0) is the constant term, coefficients(2) is the x^2 - 1 term, etc.
			 * Use the second order recurrence relation:
			 * https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation
			 * Code should be 15-20 lines
			 */
		}

		arr evaluate_hermite_polynomial_standardized(const arr& coefficients, const arr& x, const double mu, const double sigma) {
			/* Evaluate a polynomial expressed in the Hermite basis at an array of points,
			 * after standardizing them with the given mu and sigma
			 * Code should be 1-2 lines
			 */
		}

		arr2 vandermonde(const arr& x, const int degree) {
			/* Generate Vandermonde matrix at an array of points
			 * Returns 2d array of shape (x.size(), degree+1)
			 * First column is all ones
			 * Code should be 5-10 lines
			 */
		}

		arr2 hermite_vandermonde(const arr& x, const int degree) {
			/* Generate Hermite "Vandermonde" matrix at an array of points
			 * Returns 2d array of shape (x.size(), degree+1)
			 * First column is all ones
			 * Code should be 5-10 lines
			 */
		}

		std::tuple<arr2, double, double> hermite_vandermonde_standardized(const arr& x, const int degree) {
			/* Generate Hermite "Vandermonde" matrix at an array of points,
			 * standardizing them first
			 * Returns a 3-tuple containing
			 *  a 2d array of shape (x.size, degree+1)
			 *  mu
			 *  sigma
			 * Should be 5-10 lines
			 */
		}

		std::pair<arr, arr> fit_linear_regression(const arr2& X, const arr& y, const bool precondition) {
			/* Fit linear regression to solve the equation X'X beta  = X'y
			 * If precondition is set, divide each column by its l2 norm for numerical stability
			 * Returns a pair containing
			 *  coefficients beta
			 *  fitted value ~= X*beta
			 * Should be 15-20 lines
			 */
		}
	}

}

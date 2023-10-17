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

			/* We can also use Vandermonte matrix to calculate
			* arr2 mat = vandermonde(x, coefficients.size() - 1);
			arr result = arr::Zero(x.size());
			for (int i = 0; i < x.size(); ++i)
			{
				result(i) = coefficients.matrix().dot(mat.row(i).matrix());
			}
			return result;
			*/

			 // Initialize result array with the highest order coefficient
			arr result = arr::Constant(x.size(), coefficients[coefficients.size() - 1]);

			// Horner's rule
			for (int i = coefficients.size() - 2; i >= 0; --i) {
				result = result * x + coefficients[i];
			}

			return result;
		}

		arr evaluate_hermite_polynomial(const arr& coefficients, const arr& x) {
			/* Evaluate a polynomial expressed in the Hermite basis at an array of points
			 *
			 * coefficients(0) is the constant term, coefficients(2) is the x^2 - 1 term, etc.
			 * Use the second order recurrence relation:
			 * https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation
			 * Code should be 15-20 lines
			 */
			 // Initialize array to store values of the Hermite polynomials
			arr2 H = hermite_vandermonde(x, coefficients.size() - 1);
			arr result = arr::Zero(x.size());
			for (int i = 0; i < x.size(); ++i)
			{
				result(i) = coefficients.matrix().dot(H.row(i).matrix());
			}
			return result;
		}

		arr evaluate_hermite_polynomial_standardized(const arr& coefficients, const arr& x, const double mu, const double sigma) {
			/* Evaluate a polynomial expressed in the Hermite basis at an array of points,
			 * after standardizing them with the given mu and sigma
			 * Code should be 1-2 lines
			 */
			arr xhat = (x - mu) / sigma;
			return evaluate_hermite_polynomial(coefficients, xhat);
		}
		
		arr2 vandermonde(const arr& x, const uint degree) {
			/* Generate Vandermonde matrix at an array of points
			 * Returns 2d array of shape (x.size(), degree+1)
			 * First column is all ones
			 * Code should be 5-10 lines
			 */
			arr2 matrix(x.size(), degree + 1);
			for (uint i = 0; i <= degree; ++i) {
				matrix.col(i) = x.pow(i);
			}
			return matrix;
		}

        arr2 hermite_vandermonde(const arr& x, const uint degree) {
            /* Generate Hermite "Vandermonde" matrix at an array of points
             * Returns 2d array of shape (x.size(), degree+1)
             * First column is all ones
             * Code should be 5-10 lines
             */

            // use The "probabilist's Hermite polynomials"
            arr2 matrix = arr2::Zero(x.size(), degree + 1);

            // First column is all ones
            matrix.col(0) = arr::Ones(x.size());

            if (degree > 0) {
                // Second column is x
                matrix.col(1) = x;

                // Use recurrence relation for subsequent columns
                for (uint i = 2; i <= degree; ++i) {
                    matrix.col(i) = (matrix.col(i - 1) * x) - (matrix.col(i - 2) * (i - 1));
                }
            }

            return matrix;
        }

		std::tuple<arr2, double, double> hermite_vandermonde_standardized(const arr& x, const uint degree) {
			/* Generate Hermite "Vandermonde" matrix at an array of points,
			 * standardizing them first
			 * Returns a 3-tuple containing
			 *  a 2d array of shape (x.size, degree+1)
			 *  mu
			 *  sigma
			 * Should be 5-10 lines
			 */
			double mean = x.mean();
			double std_dev = std::sqrt((x - mean).square().sum() / x.size());

			// Standardize x
			arr x_standardized = (x - mean) / std_dev;

			// Generate Hermite Vandermonde matrix with standardized x
			arr2 matrix = hermite_vandermonde(x_standardized, degree);

			return { matrix, mean, std_dev };
		}

		std::pair<arr,arr> fit_linear_regression(const arr2& X, const arr& y, const bool precondition) {
			/* Fit linear regression to solve the equation X'X beta  = X'y
			 * If precondition is set, divide each column by its l2 norm for numerical stability
			 * Returns a pair containing
			 *  coefficients beta
			 *  fitted value ~= X*beta
			 * Should be 15-20 lines
			 */
			arr2 D = ArrayXXd::Zero(X.cols(), X.cols());
			if (precondition) {
				for (int i = 0; i < X.cols(); ++i) {
					double norm = 0.0;
					for (int j = 0; j < X.rows(); ++j) {
						norm += X(j, i) * X(j, i);
					}
					norm = std::sqrt(norm);

					if (norm > 0) {
						D(i, i) = norm;
					}
				}
			}
			Eigen::MatrixXd X_precond = X.matrix() * D.matrix().inverse();
			Eigen::BDCSVD<Eigen::MatrixXd> X_svd = X_precond.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
			Eigen::MatrixXd V = X_svd.matrixV();
			Eigen::MatrixXd U = X_svd.matrixU();
			Eigen::MatrixXd sigma = X_svd.singularValues().asDiagonal();
			Eigen::MatrixXd part = (U.transpose()) * (y.matrix());
			arr beta1 = D.matrix().inverse() * V * sigma.inverse() * part;
			arr fitted = U * part;

			return std::make_pair(beta1, fitted);
		}
	}

}

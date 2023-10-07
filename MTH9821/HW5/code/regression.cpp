#include "regression.h"
#include <tuple>

namespace linalg {

	namespace regression {
		
		using namespace Eigen;

		arr evaluate_polynomial(const arr& coefficients, const arr& x) {
    			unsigned int len = x.size();
    			arr result(len);
                for (int i = 0; i < x.size(); ++i) {
                    double sum = 0;
                    for (int j = coefficients.size() - 1; j >= 0; --j) {
                        sum = sum * x(i) + coefficients(j);
                    }
                    result(i) = sum;
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
            //Chat GPT
            int n = coefficients.size();
            Eigen::ArrayXd result = Eigen::ArrayXd::Zero(x.size());
            
            if(n == 0) return result;
            
            for(int i = 0; i < x.size(); ++i) {
                double h_prev = 1.0;
                double h_curr = 2 * x(i) - 1;
                double temp;

                if(n > 1) {
                    for(int j = 2; j < n; ++j) {
                        temp = h_curr;
                        h_curr = 2 * (x(i) * h_curr - (j - 1) * h_prev);
                        h_prev = temp;
                    }
                }

                result(i) = coefficients(0) * h_prev + (n > 1 ? coefficients(1) * h_curr : 0);

                for(int j = 2; j < n; ++j) {
                    result(i) += coefficients(j) * h_curr;
                }
            }
            
            return result;
		}

		arr evaluate_hermite_polynomial_standardized(const arr& coefficients, const arr& x, const double mu, const double sigma) {
			/* Evaluate a polynomial expressed in the Hermite basis at an array of points,
			 * after standardizing them with the given mu and sigma
			 * Code should be 1-2 lines
			 */
            return evaluate_hermite_polynomial(coefficients, (x - mu) / sigma);
		}
		
		arr2 vandermonde(const arr& x, const uint degree) {
			/* Generate Vandermonde matrix at an array of points
			 * Returns 2d array of shape (x.size(), degree+1)
			 * First column is all ones
			 * Code should be 5-10 lines
			 */
            int n = x.size();
            MatrixXd V = MatrixXd::Ones(n, degree + 1); // create matrix with ones
            for (int i = 0; i < n; ++i) {
                for (unsigned int j = 1; j <= degree; ++j) {
                    V(i, j) = std::pow(x(i), j); // compute powers
                }
            }
            return V;
		}

		arr2 hermite_vandermonde(const arr& x, const uint degree) {
            arr2 result(x.size(), degree + 1);
            result.col(0) = arr::Ones(x.size()); // First column is all ones

            for (uint j = 1; j <= degree; ++j) {
                arr coefficients = arr::Zero(j + 1);
                coefficients(j) = 1; // Set coefficient for x^j term to 1, others are 0
                result.col(j) = evaluate_hermite_polynomial(coefficients, x);
            }

            return result;
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
            double mu = x.mean();
            double sigma = std::sqrt((x - mu).square().mean());

            arr standardized_x = (x - mu) / sigma;
            arr2 result = hermite_vandermonde(standardized_x, degree);

            return std::make_tuple(result, mu, sigma);
		}

		std::pair<arr,arr> fit_linear_regression(const arr2& X, const arr& y, const bool precondition) {
			/* Fit linear regression to solve the equation X'X beta  = X'y
			 * If precondition is set, divide each column by its l2 norm for numerical stability
			 * Returns a pair containing
			 *  coefficients beta
			 *  fitted value ~= X*beta
			 * Should be 15-20 lines
			 */
            arr2 X_copy = X;
            if(precondition) {
                for(int i = 0; i < X.cols(); ++i) {
                    double norm = X_copy.col(i).matrix().norm();
                    if(norm != 0) { // avoid division by zero
                        X_copy.col(i) /= norm;
                    }
                }
            }
            
            // Solving the linear regression problem using JacobiSVD for numerical stability
            Eigen::JacobiSVD<arr2> svd(X_copy, Eigen::ComputeThinU | Eigen::ComputeThinV);
            arr beta = svd.solve(y);

            // If preconditioned, adjust beta accordingly
            if(precondition) {
                for(int i = 0; i < X.cols(); ++i) {
                    beta(i) /= X.col(i).norm();
                }
            }

            arr fitted = X * beta;
            return std::make_pair(beta, fitted);
		}
	}

}

#include "regression.h"
#include <tuple>
#include <utility>  // for std::pair
//#include <Eigen/SVD>



namespace linalg {

	namespace regression {
		
		using namespace Eigen;

		Eigen::VectorXd evaluate_polynomial(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& x) {
                //Test case is written
    			unsigned int len = x.size();
    			Eigen::VectorXd result(len);
                for (int i = 0; i < x.size(); ++i) {
                    double sum = 0;
                    for (int j = coefficients.size() - 1; j >= 0; --j) {
                        sum = sum * x(i) + coefficients(j);
                    }
                    result(i) = sum;
                }
    			return result;
		}

		Eigen::VectorXd evaluate_hermite_polynomial(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& x) {
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

		Eigen::VectorXd evaluate_hermite_polynomial_standardized(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& x, const double mu, const double sigma) {
			/* Evaluate a polynomial expressed in the Hermite basis at an array of points,
			 * after standardizing them with the given mu and sigma
			 * Code should be 1-2 lines
			 */
            Eigen::VectorXd standardized_x = (x.array() - mu) / sigma;
            return evaluate_hermite_polynomial(coefficients, standardized_x);
		}
		
		Eigen::MatrixXd vandermonde(const Eigen::VectorXd& x, const uint degree) {
            //Vandermode tested worked
			/* Generate Vandermonde matrix at an array of points
			 * Returns 2d array of shape (x.size(), degree+1)
			 * First column is all ones
			 * Code should be 5-10 lines
			 */
            int n = x.size();
            Eigen::MatrixXd V = Eigen::MatrixXd::Ones(n, degree + 1); // create matrix with ones
            for (int i = 0; i < n; ++i) {
                for (unsigned int j = 1; j <= degree; ++j) {
                    V(i, j) = std::pow(x(i), j); // compute powers
                }
            }
            return V;
		}

		Eigen::MatrixXd hermite_vandermonde(const Eigen::VectorXd& x, const uint degree) {
            Eigen::MatrixXd result(x.size(), degree + 1);
            result.col(0) = Eigen::VectorXd::Ones(x.size()); // First column is all ones

            for (uint j = 1; j <= degree; ++j) {
                Eigen::VectorXd coefficients = arr::Zero(j + 1);
                coefficients(j) = 1; // Set coefficient for x^j term to 1, others are 0
                result.col(j) = evaluate_hermite_polynomial(coefficients, x);
            }

            return result;
		}

		std::tuple<Eigen::MatrixXd, double, double> hermite_vandermonde_standardized(const Eigen::VectorXd& x, const uint degree) {
			/* Generate Hermite "Vandermonde" matrix at an array of points,
			 * standardizing them first
			 * Returns a 3-tuple containing
			 *  a 2d array of shape (x.size, degree+1)
			 *  mu
			 *  sigma
			 * Should be 5-10 lines
			 */
            double mu = x.mean();

            double sigma = std::sqrt((x.array() - mu).square().mean());

            Eigen::VectorXd standardized_x = (x.array() - mu) / sigma;
            Eigen::MatrixXd result = hermite_vandermonde(standardized_x, degree);

            return std::make_tuple(result, mu, sigma);
		}

		std::pair<Eigen::VectorXd, Eigen::VectorXd> fit_linear_regression(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const bool precondition) {
            Eigen::MatrixXd X_modified = X; // Create a copy of X to modify if precondition is set
            
            if (precondition) {
                // Precondition X: Divide each column by its L2 norm
                for (int i = 0; i < X.cols(); ++i) {
                    double norm = X.col(i).matrix().norm();
                    if (norm > 0) {
                        X_modified.col(i) /= norm;
                    }
                }
            }
            
             // Calculate beta using the normal equation without explicitly inverting
            Eigen::MatrixXd Xty = X_modified.transpose() * y;
            Eigen::MatrixXd XtX = X_modified.transpose() * X_modified;

            // Solve the linear system instead of inverting
            Eigen::VectorXd beta = XtX.colPivHouseholderQr().solve(Xty);

            // Calculate the fitted values
            Eigen::VectorXd fitted_values = X * beta.matrix();

            return std::make_pair(beta, fitted_values);
        }
    }
}

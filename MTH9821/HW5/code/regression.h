#ifndef REGRESSION_H
#define REGRESSION_H

#include <tuple>
#include <Eigen/Dense>

namespace linalg {
    namespace regression {
		typedef Eigen::ArrayXd arr;
		typedef Eigen::ArrayXXd mat;
        Eigen::VectorXd evaluate_polynomial(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& x);
        Eigen::VectorXd evaluate_hermite_polynomial(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& x);
        Eigen::VectorXd evaluate_hermite_polynomial_standardized(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& z, const double mu, const double sigma);
        Eigen::MatrixXd vandermonde(const Eigen::VectorXd& x, const uint degree);
        Eigen::MatrixXd hermite_vandermonde(const Eigen::VectorXd& x, const uint degree);
        std::tuple<Eigen::MatrixXd, double, double> hermite_vandermonde_standardized(const Eigen::VectorXd& x, const uint degree);
        std::pair<Eigen::VectorXd,Eigen::VectorXd> fit_linear_regression(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const bool precondition) ;
    }
}

#endif // REGRESSION_H

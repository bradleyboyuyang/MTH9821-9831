#ifndef MC_REGRESSION_H
#define MC_REGRESSION_H

#include <vector>
#include <Eigen/Dense>

namespace montecarlo {
    namespace regression {
		typedef Eigen::ArrayXd arr;
		typedef Eigen::ArrayXXd arr2;

        struct MCRegression {
            const uint _times;
            double _coeff;
            MCRegression(const uint times);
            // for t = 0
            double fit_predict_at_0(const arr& y);
            double predict_at_0() const;
            // for t > 0
            virtual arr fit_predict(const uint time, const arr& x, const arr& y) = 0;
            virtual arr predict(const uint time, const arr& x) const = 0;
            virtual ~MCRegression();
        };

        struct PolynomialMCRegression: MCRegression {
            const uint degree;
            const bool precondition;
            std::vector<arr> _coeffs;
            PolynomialMCRegression(const uint times, const uint degree, const bool precondition);
            arr fit_predict(const uint time, const arr& x, const arr& y) override;
            arr predict(const uint time, const arr& x) const override;
        };

        struct HermiteMCRegression: MCRegression {
            const uint degree;
            const bool precondition;
            const bool standardize;
            std::vector<arr> _coeffs;
            std::vector<double> _means;
            std::vector<double> _stds;
            HermiteMCRegression(const uint times, const uint degree, const bool precondition, const bool standardize);
            arr fit_predict(const uint time, const arr& x, const arr& y) override;
            arr predict(const uint time, const arr& x) const override;
        };
    }
}

#endif // MC_REGRESSION_H

//
//  finite_difference_pricer.hpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/11/22.
//

#ifndef finite_difference_pricer_hpp
#define finite_difference_pricer_hpp

#include "EuropeanOption.hpp"
#include <vector>

enum Euler {
    Explicit,
    Implicit,
    ImEx,
};

enum OptionType {
    European,
    American,
};

class FiniteDifferencePricer {
protected:
    // Option data
    double S0_;
    double K_;
    double T_;
    double sigma_;
    double r_;
    double q_;
    
    // Heat equation transformation coefficients
    double a_;
    double b_;
    
    // Finite difference parameters
    double tau_final_;
    double x_l_;
    double x_r_;
//    double dt_;
//    double dx_;
//    std::size_t N_;
//    double alpha_;
    
    // Boundary conditions
    std::function<double (double)> boundary_tau_0_;
    std::function<double (double)> boundary_x_l_;
    std::function<double (double)> boundary_x_r_;
    
    // Modules used in finite difference pricing
    // 1. Computational domain
    std::tuple<std::size_t, double, double, double> DomainParams(std::size_t M) const;
    
    std::pair<std::vector<double>, std::size_t> BuildMesh(std::size_t N, double dx) const;
    
    // 2. Boundary conditions
    void GenerateBoundaryConditions(const OptionType& option_type);
    
    // 3. Finite difference scheme
    std::vector<std::vector<double>> FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau, const Euler& euler, const OptionType& option_type) const;
    
    // 4. Pointwise convergence
    std::vector<double> FindValue(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const;
    
    // 5. RMS error
    double FindRmsError(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) const;
    
    // 6. Greeks
    std::vector<double> FindGreeks(std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const;
    
    // 7. Variance reduction
    double VarianceReduction(double V_approx, std::size_t M, const Euler& euler) const;
    
public:
    FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q);
    FiniteDifferencePricer(const EuropeanOption& option);
    ~FiniteDifferencePricer() = default;
    void PrintVector(const std::vector<double>& vec) const;
    
    // Finite difference hyperparameters
    static double alpha_temp_;
    static std::size_t M_init_;
    
    std::vector<double> Price(std::size_t M, const Euler& euler, const OptionType& option_type);
    std::vector<double> Price(std::size_t M, std::size_t N, const Euler& euler, const OptionType& option_type);
    
    friend class FiniteDifferenceEngine;
    
};


#endif /* finite_difference_pricer_hpp */

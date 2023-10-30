//
//  finite_difference_pricer_discrete_dividend.hpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/18/22.
//

#ifndef finite_difference_pricer_discrete_dividend_hpp
#define finite_difference_pricer_discrete_dividend_hpp

#include "finite_difference_pricer.hpp"
#include "EigenCommonHeader.h"

// ERSATZ. Only handles one discrete dividend event.

class FiniteDifferencePricerDiscreteDividend {
private:
    // Option data
    double S0_;
    double K_;
    double T_;
    double sigma_;
    double r_;
    double q_;

    // Proportional dividend (discrete)
    double prop_;
    
    // Heat equation transformation coefficients
    double a_;
    double b_;
    
    // Time of dividend (just once)
    // Change to vector if there are multiple dividends
    double t_div_;
    double tau_div_;
    
    double tau_final_;
    
    // Boundary x vals
    double x_l_;
    double x_r_;
    double x_compute_;
    
    // Number of nodes per level
    std::vector<std::size_t> N_;
    
    // Boundary conditions
    std::function<double (double)> boundary_tau_0_;
    std::function<double (double)> boundary_x_l_;
    std::function<double (double)> boundary_x_r_;
    
public:
    FiniteDifferencePricerDiscreteDividend(double S0, double K, double T, double sigma, double r, double q, double t_div, double prop);
    FiniteDifferencePricerDiscreteDividend(const EuropeanOption& option, double t_div, double prop);
    ~FiniteDifferencePricerDiscreteDividend() = default;
    void PrintVector(const std::vector<double>& vec) const;
    
    // Finite difference hyperparameters
    static double alpha_temp_;
    static std::size_t M_init_;
    
    std::vector<double> PriceCall(std::size_t M, const Euler& euler, const OptionType& option_type);
    
private:
    // Modules used in finite difference pricing
    // 1. Computational domain
    std::pair<std::vector<double>, std::size_t> BuildMesh(std::size_t M_1);
    
    // 2. Boundary conditions
    void GenerateBoundaryConditions(const OptionType& option_type);
    
    // 3. Finite difference scheme
    std::vector<std::vector<double>> FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::vector<double> u_mesh, std::size_t M, double dtau, const Euler& euler, const OptionType& option_type, double curr_tau, double tau_end) const;
    
//    // 4. Pointwise convergence
//    std::vector<double> FindValue(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const;
    
    // 6. Greeks
    std::vector<double> FindGreeks(std::size_t x_node, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const;
    
    std::vector<std::vector<double>> FullExplAmer(double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, std::size_t M, double dtau, double curr_tau, double tau_end) const;
    std::vector<std::vector<double>> FullImExAmer(double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, std::size_t M, double dtau, double curr_tau, double tau_end) const;
    
    void StepExplAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    
    void StepImExAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& A, const mat& b_multiplier) const;
};


#endif /* finite_difference_pricer_discrete_dividend_hpp */

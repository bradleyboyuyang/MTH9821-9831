//
//  finite_difference_pricer_discrete_dividend.cpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/18/22.
//

#include "finite_difference_pricer_discrete_dividend.hpp"
#include <cmath>
#include <iostream>
#include <exception>

#include "LinearSolver.hpp"

double FiniteDifferencePricerDiscreteDividend::alpha_temp_ = 0.45;
std::size_t FiniteDifferencePricerDiscreteDividend::M_init_ = 4;

void FiniteDifferencePricerDiscreteDividend::PrintVector(const std::vector<double>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << '\t';
    }
    std::cout << std::endl;
}

FiniteDifferencePricerDiscreteDividend::FiniteDifferencePricerDiscreteDividend(double S0, double K, double T, double sigma, double r, double q, double t_div, double prop) : t_div_(t_div), prop_(prop), S0_(S0), K_(K), T_(T), sigma_(sigma), q_(q), r_(r) {
    
    assert(q == 0);
    
    //  1. Computational domain (part 1)
    double sigma2 = sigma * sigma;

    tau_final_ = T * sigma2 / 2.;
    tau_div_ = (T - t_div) * sigma2 / 2.;

    // temp x_l1_ and x_r1_
    x_l_ = std::log(S0 / K) + (r - q - sigma2 / 2.) * T - 3. * sigma * std::sqrt(T);
    x_r_ = std::log(S0 / K) + (r - q - sigma2 / 2.) * T + 3. * sigma * std::sqrt(T);

    x_compute_ = log(S0 / K) + log(1. - prop);
    
    // Heat equation transformation coefficients
    a_ = (r - q) / sigma2 - .5;
    b_ = (a_ + 1.) * (a_ + 1.) + 2. * q / sigma2;
}

FiniteDifferencePricerDiscreteDividend::FiniteDifferencePricerDiscreteDividend(const EuropeanOption& option, double t_div, double prop) : FiniteDifferencePricerDiscreteDividend(option.S_, option.K_, option.T_, option.sigma_, option.r_, option.q_, t_div, prop) {}

// MARK: Pricing function exposed to user
std::vector<double> FiniteDifferencePricerDiscreteDividend::PriceCall(std::size_t M, const Euler& euler, const OptionType& option_type) {
    
    std::vector<double> res;
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t x_node; // Node where x0 is contained
    std::tie(x_mesh, x_node) = this->BuildMesh(M);
    double alpha = FiniteDifferencePricerDiscreteDividend::alpha_temp_;
    
    // Here, x_l and x_r are modified!!!
    this->GenerateBoundaryConditions(option_type);
    
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    
    double dtau_1 = tau_div_ / M;
    double dx = std::sqrt(dtau_1 / alpha);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference(alpha, x_mesh, u_mesh, M, dtau_1, euler, option_type, 0., tau_div_);
    u_mesh = u_meshes[0];
    
    auto AdjustXDividend = [=](double x)->double {
        return x - std::log(1. - prop_);
    };
    std::transform(x_mesh.cbegin(), x_mesh.cend(), x_mesh.begin(), AdjustXDividend);
//    std::cout << x_l_ << '\t' << x_r_ << std::endl;
    x_l_ = AdjustXDividend(x_l_);
    x_r_ = AdjustXDividend(x_r_);
//    std::cout << x_l_ << '\t' << x_r_ << std::endl;
    
    auto AdjustUDividend = [=](double u)->double {
        return std::pow(1. - prop_, -a_) * u;
    };
    std::transform(u_mesh.cbegin(), u_mesh.cend(), u_mesh.begin(), AdjustUDividend);
    
    std::size_t M_2 = std::ceil((tau_final_ - tau_div_) / (alpha * dx * dx));
    double dtau_2 = (tau_final_ - tau_div_) / double(M_2);
    double alpha_2 = dtau_2 / (dx * dx);
    
//    std::cout << M_2 << '\t' << alpha_2 << std::endl;
//    std::cout << dtau_1 << '\t' << dtau_2 << std::endl;
//    std::cout << tau_div_ << tau_final_ << std::endl;
    
    u_meshes = this->FiniteDifference(alpha_2, x_mesh, u_mesh, M_2, dtau_2, euler, option_type, tau_div_, tau_final_);
    u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    res.push_back(u_mesh[x_node]);
    
    // 4. Pointwise convergence
    double V_approx = std::exp(-a_ * x_mesh[x_node] - b_ * tau_final_) * u_mesh[x_node];
    res.push_back(V_approx);
    
    // 6. Greeks
    std::vector<double> Greeks = this->FindGreeks(x_node, x_mesh, u_mesh, u_mesh_prev, dtau_2, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::pair<std::vector<double>, std::size_t> FiniteDifferencePricerDiscreteDividend::BuildMesh(std::size_t M_1) {
    // More finite difference parameters
    double dtau_1 = tau_div_ / M_1;
    // N: number of x intervals on the x-axis
    double dx = std::sqrt(dtau_1 / FiniteDifferencePricerDiscreteDividend::alpha_temp_);
    
    std::size_t N_l = std::ceil((x_compute_ -  x_l_) / dx);
    std::size_t N_r = std::ceil((x_r_ -  x_compute_) / dx);
    
    x_l_ = x_compute_ - N_l * dx;
    x_r_ = x_compute_ + N_r * dx;
    
    std::vector<double> x_mesh(N_l + N_r + 1);
    for (double i = 0; i < x_mesh.size(); i++) {
        x_mesh[i] = x_compute_ + dx * (i - N_l);
    }
    
    assert(x_mesh.front() == x_l_);
    assert(x_mesh.back() == x_r_);

    return std::make_pair(x_mesh, N_l);
}

void FiniteDifferencePricerDiscreteDividend::GenerateBoundaryConditions(const OptionType& option_type) {
    switch (option_type) {
//        case European:
//            // !!!: Not modified.
//
//            // 2. Boundary conditions
//            boundary_tau_0_ = [=](double x)->double {
//                // K * exp(ax) * max(1 - exp(x), 0)
//                if (x < 0.) {
//                    return K_ * std::exp(a_ * x) * (1. - std::exp(x));
//                } else {
//                    return 0.;
//                }
//            };
//
//            boundary_x_l_ = [=](double tau)->double {
//                // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
//                return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
//            };
//
//            boundary_x_r_ = [](double tau)->double {
//                return 0.;
//            };
//            break;
        default:
            // 2. Boundary conditions
            boundary_tau_0_ = [=](double x)->double {
                // K * exp(ax) * max(1 - exp(x), 0)
//                std::cout << "Called" << std::endl;
                if (x > 0.) {
                    // MARK: Boundary condition for American call!!
                    return K_ * std::exp(a_ * x) * (std::exp(x) - 1.);
                } else {
                    return 0.;
                }
            };
            
            boundary_x_l_ = [=](double tau)->double {
//                return K_ * std::exp(a_ * x_l_ + b_ * tau) * (1. - std::exp(x_l_));
                if (tau <= tau_div_) {
                    return 0.;
                } else {
                    return 0.;
                }
            };
            
            boundary_x_r_ = [=](double tau)->double {
//                return 0.;
                if (tau <= tau_div_) {
                    return K_ * std::exp(a_ * x_r_ + b_ * tau) * (std::exp(x_r_) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
                } else {
//                    return K_ * std::exp(a_ * x_r_ + b_ * tau) * (std::exp(x_r_) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
//                    std::cout << "HERE" << std::endl;
//                    return K_ * std::exp(a_ * (x_r_ - std::log(1. - prop_)) + b_ * tau) * (std::exp(x_r_) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
                    return K_ * std::exp(a_ * (x_r_) + b_ * tau) * (std::exp(x_r_ + std::log(1. - prop_)) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
                }
                
            };
            break;
    }
}

std::vector<std::vector<double>> FiniteDifferencePricerDiscreteDividend::FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::vector<double> u_mesh, std::size_t M, double dtau, const Euler& euler, const OptionType& option_type, double curr_tau, double tau_end) const {
    assert(option_type == American);
    
    switch (euler) {
        case Explicit:
            return this->FullExplAmer(alpha, x_mesh, u_mesh, M, dtau, curr_tau, tau_end);
            break;
        case ImEx:
            return this->FullImExAmer(alpha, x_mesh, u_mesh, M, dtau, curr_tau, tau_end);
        default:
            return this->FullExplAmer(alpha, x_mesh, u_mesh, M, dtau, curr_tau, tau_end);
            break;
    }
}

std::vector<std::vector<double>> FiniteDifferencePricerDiscreteDividend::FullExplAmer(double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, std::size_t M, double dtau, double curr_tau, double tau_end) const {
//    this->PrintVector(u_mesh);
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        curr_tau += dtau;
        this->StepExplAmer(curr_tau, alpha, x_mesh, u_mesh);
//        this->PrintVector(u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepExplAmer(tau_end, alpha, x_mesh, u_mesh);
//    this->PrintVector(u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferencePricerDiscreteDividend::FullImExAmer(double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, std::size_t M, double dtau, double curr_tau, double tau_end) const {
    
    mat A(mat::Zero(x_mesh.size() - 2, x_mesh.size() - 2));   // Initialize with zero matrix
    A(0, 0) = 1. + alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha * .5;
        A(i, i - 1) = -alpha * .5;
        A(i, i) = 1. + alpha;
    }
    
    // Build matrix needed for the construction of b
    mat b_multiplier(mat::Zero(A.rows(), A.cols()));   // Initialize with zero matrix
    
    b_multiplier(0, 0) = 1. - alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        b_multiplier(i - 1, i) = alpha * .5;
        b_multiplier(i, i - 1) = alpha * .5;
        b_multiplier(i, i) = 1. - alpha;
    }
//    this->PrintVector(u_mesh);
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        curr_tau += dtau;
        this->StepImExAmer(curr_tau, alpha, x_mesh, u_mesh, A, b_multiplier);
//        this->PrintVector(u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepImExAmer(tau_end, alpha, x_mesh, u_mesh, A, b_multiplier);
//    this->PrintVector(u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

void FiniteDifferencePricerDiscreteDividend::StepExplAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(boundary_x_l_(tau));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        // Get the corresponding European option's value
        double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
        
        // Find early exercise
        double early_ex_premium = 0.;
        
        if (x_mesh[pos] > 0.) {
            early_ex_premium = K_ * std::exp(a_ * x_mesh[pos]) * (std::exp(x_mesh[pos]) - 1.);
        }
        
        // Compare and add to mesh
        
        new_u_mesh.push_back(std::max(euro_val, early_ex_premium));
    }
    
    // Right boundary
    new_u_mesh.push_back(boundary_x_r_(tau));
    
    u_mesh = std::move(new_u_mesh);
}

void FiniteDifferencePricerDiscreteDividend::StepImExAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& A, const mat& b_multiplier) const {
    
    // Prepare b
    vec b(u_mesh.size() - 2);
    
    // [Some matrix] * u
    vec u(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, u.begin());
    
    b = b_multiplier * u;
    
    // Boundary conditions
    b(0) += boundary_x_l_(tau) * alpha * .5;
    b(b.size() - 1) += boundary_x_r_(tau) * alpha * .5;
    b(0) += u_mesh.front() * alpha * .5;
    b(b.size() - 1) += u_mesh.back() * alpha * .5;
    
//    double old_right_boundary = *(u_mesh.rbegin());
    
    // Add boundary condition
    u_mesh.front() = boundary_x_l_(tau);
    u_mesh.back() = boundary_x_r_(tau);
    
    auto find_early_ex = [=](double x)->double {
        if (x > 0.) {
            return K_ * std::exp(a_ * x) * (std::exp(x) - 1.);
        } else {
            return 0.;
        }
    };
    
    vec early_ex_premium(x_mesh.size() - 2);
    std::transform(x_mesh.cbegin() + 1, x_mesh.cend() - 1, early_ex_premium.begin(), find_early_ex);

    // Solve linear system
    LinearSolver solver;
    vec sol = solver.LUSolve(A, b);
    
    std::copy(sol.cbegin(), sol.cend(), u_mesh.begin() + 1);
}

std::vector<double> FiniteDifferencePricerDiscreteDividend::FindGreeks(std::size_t x_node, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const {
    
    std::vector<double> greeks;
    
    double x_small = x_mesh[x_node - 1];
    double x_mid = x_mesh[x_node];
    double x_large = x_mesh[x_node + 1];
    
    double s_small = K_ * std::exp(x_small);
    double s_mid = K_ * std::exp(x_mid);
    double s_large = K_ * std::exp(x_large);
    
    double v_small = std::exp(-a_ * x_small - b_ * tau_final_) * u_mesh[x_node - 1];
    double v_mid = std::exp(-a_ * x_mid - b_ * tau_final_) * u_mesh[x_node];
    double v_large = std::exp(-a_ * x_large - b_ * tau_final_) * u_mesh[x_node + 1];
    
    double delta = (v_large - v_small) / (s_large - s_small);
    double gamma = ((s_mid - s_small) * v_large - (s_large - s_small) * v_mid + (s_large - s_mid) * v_small) / ((s_mid - s_small) * (s_large - s_mid) * (s_large - s_small) * .5);
    
    greeks.push_back(delta);
    greeks.push_back(gamma);
    
    double dt = 2. * dtau / (sigma_ * sigma_);
    double v_prev = std::exp(-a_ * x_mid - b_ * (tau_final_ - dtau)) * u_mesh_prev[x_node];
    double theta = (v_prev - v_mid) / dt;
    greeks.push_back(theta);
    
    return greeks;
}

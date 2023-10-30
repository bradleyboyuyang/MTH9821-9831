//
//  finite_difference_pricer.cpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/11/22.
//

#include "finite_difference_pricer.hpp"
#include <cmath>
#include <iostream>
#include <exception>

#include "finite_difference_engine.hpp"

double FiniteDifferencePricer::alpha_temp_ = 0.45;
std::size_t FiniteDifferencePricer::M_init_ = 4;

FiniteDifferencePricer::FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q) : S0_(S0), K_(K), T_(T), sigma_(sigma), q_(q), r_(r) {
    
    //  1. Computational domain (part 1)
    double sigma2 = sigma * sigma;
    
    tau_final_ = T * sigma2 / 2.;
    
    x_l_ = std::log(S0 / K) + (r - q - sigma2 / 2.) * T - 3. * sigma * std::sqrt(T);
    x_r_ = std::log(S0 / K) + (r - q - sigma2 / 2.) * T + 3. * sigma * std::sqrt(T);
    
    // Heat equation transformation coefficients
    a_ = (r - q) / sigma2 - .5;
    b_ = (a_ + 1.) * (a_ + 1.) + 2. * q / sigma2;
}

FiniteDifferencePricer::FiniteDifferencePricer(const EuropeanOption& option) : FiniteDifferencePricer(option.S_, option.K_, option.T_, option.sigma_, option.r_, option.q_) {}

void FiniteDifferencePricer::PrintVector(const std::vector<double>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << '\t';
    }
    std::cout << std::endl;
}

// MARK: Pricing function exposed to user

std::vector<double> FiniteDifferencePricer::Price(std::size_t M, std::size_t N, const Euler& euler, const OptionType& option_type) {
    
    std::vector<double> res;
    
    double dtau = tau_final_ / M;
    double dx = (x_r_ - x_l_) / N;
    double alpha = dtau / (dx * dx);
//    std::cout << alpha << std::endl;
    
    this->GenerateBoundaryConditions(option_type);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i; // Interval where x_mesh is contained
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    // Finite difference scheme
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference(alpha, x_mesh, M, dtau, euler, option_type);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->FindValue(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    return res;
    
}

std::vector<double> FiniteDifferencePricer::Price(std::size_t M, const Euler& euler, const OptionType& option_type) {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    this->GenerateBoundaryConditions(option_type);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i; // Interval where x_mesh is contained
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    // Finite difference scheme
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference(alpha, x_mesh, M, dtau, euler, option_type);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->FindValue(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 5. RMS error
    if (option_type == European) {
        double error_RMS = this->FindRmsError(x_mesh, u_mesh);
        res.push_back(error_RMS);
    }
    
    // 6. Greeks
    std::vector<double> Greeks = this->FindGreeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    // 7. Variance reduction
    if (option_type == American) {
        double variance_reduced_value = this->VarianceReduction(V_approx, M, euler);
        res.push_back(variance_reduced_value);
    }
    
    return res;
}

// MARK: Modules used in finite difference pricing

std::tuple<std::size_t, double, double, double> FiniteDifferencePricer::DomainParams(std::size_t M) const {
    
    // More finite difference parameters
    double dtau = tau_final_ / M;
    // N: number of x intervals on the x-axis
    std::size_t N = std::floor((x_r_ - x_l_) / std::sqrt(dtau / alpha_temp_));
//    std::cout << N << std::endl;
    double dx = (x_r_ - x_l_) / N;
    double alpha = dtau / (dx * dx);
    
    return std::make_tuple(N, dtau, dx, alpha);
    
}

void FiniteDifferencePricer::GenerateBoundaryConditions(const OptionType& option_type) {
    switch (option_type) {
        case European:
            // 2. Boundary conditions
            boundary_tau_0_ = [=](double x)->double {
                // K * exp(ax) * max(1 - exp(x), 0)
                if (x < 0.) {
                    return K_ * std::exp(a_ * x) * (1. - std::exp(x));
                } else {
                    return 0.;
                }
            };
            
            boundary_x_l_ = [=](double tau)->double {
                // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
                return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
            };
            
            boundary_x_r_ = [](double tau)->double {
                return 0.;
            };
            break;
        case American:
            // 2. Boundary conditions
            boundary_tau_0_ = [=](double x)->double {
                // K * exp(ax) * max(1 - exp(x), 0)
                if (x < 0.) {
                    return K_ * std::exp(a_ * x) * (1. - std::exp(x));
                } else {
                    return 0.;
                }
            };
            
            boundary_x_l_ = [=](double tau)->double {
                // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
                return K_ * std::exp(a_ * x_l_ + b_ * tau) * (1. - std::exp(x_l_));
            };
            
            boundary_x_r_ = [](double tau)->double {
                return 0.;
            };
            break;
    }
}

std::pair<std::vector<double>, std::size_t> FiniteDifferencePricer::BuildMesh(std::size_t N, double dx) const {
    
    // Fill x mesh
    std::vector<double> x_mesh({x_l_});
    for (std::size_t i = 0; i < N; i++) {
        x_mesh.push_back(x_mesh.back() + dx);
    }
    
    // Find actual x on x_mesh
    double x_compute = std::log(S0_ / K_);
    auto x_large_it = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), x_large_it) - 1;
    
    return std::make_pair(x_mesh, interval_i);
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau, const Euler& euler, const OptionType& option_type) const {
    
    FiniteDifferenceEngine engine(this);
    return engine.RunScheme(alpha, x_mesh, M, dtau, euler, option_type);
}

std::vector<double> FiniteDifferencePricer::FindValue(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const {
    
    std::vector<double> approximations;
    
    // Find the interval containing x_compute
    double x_compute = std::log(S0_ / K_);
    auto x_large_it = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), x_large_it) - 1;
    
    // Approximation method 1
    double S_small = K_ * std::exp(x_mesh[interval_i]);
    double S_large = K_ * std::exp(x_mesh[interval_i + 1]);
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_approx = ((S_large - S0_) * V_small + (S0_ - S_small) * V_large) / (S_large - S_small);
    approximations.push_back(V_approx);
    
    // Approximation method 2
    double u_approx = ((x_mesh[interval_i + 1] - x_compute) * u_mesh[interval_i] + (x_compute - x_mesh[interval_i]) * u_mesh[interval_i + 1]) / (dx);
    double V_approx_2 = std::exp(-a_ * x_compute - b_ * tau_final_) * u_approx;
    approximations.push_back(V_approx_2);
    
    return approximations;
}

double FiniteDifferencePricer::FindRmsError(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) const {
    
    // Calculate vector of FD and BS option values
    auto V_mesh_approx_gen = [&](double x, double u)->double {
        return std::exp(-a_ * x - b_ * tau_final_) * u;
    };
    auto V_mesh_exact_gen = [&](double x)->double {
        EuropeanOption option(0., K_ * std::exp(x), K_, T_, sigma_, r_, q_);
        return option.Put();
    };
    
    std::vector<double> V_mesh_approx(x_mesh.size());
    std::vector<double> V_mesh_exact(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.cbegin(), V_mesh_approx.begin(), V_mesh_approx_gen);
    std::transform(x_mesh.cbegin(), x_mesh.cend(), V_mesh_exact.begin(), V_mesh_exact_gen);
    
    // Find RMS
    double error_sq = 0.;
    int error_count = 0;
    auto V_mesh_approx_it = V_mesh_approx.cbegin();
    auto V_mesh_exact_it = V_mesh_exact.cbegin();
    while (V_mesh_exact_it != V_mesh_exact.cend()) {
        double V_BS = *V_mesh_exact_it;
        double V_FD = *V_mesh_approx_it;
        if (V_BS > 0.00001 * S0_) {
            error_count++;
            error_sq += (V_BS - V_FD) * (V_BS - V_FD) / (V_BS * V_BS);
        }
        V_mesh_approx_it++;
        V_mesh_exact_it++;
    }
    double error_RMS = std::sqrt(error_sq / error_count);
    
    return error_RMS;
}

std::vector<double> FiniteDifferencePricer::FindGreeks(std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const {
    
    // DELTA & GAMMA
    // Get S and V of interest
    double S_smaller = K_ * std::exp(x_mesh[interval_i - 1]);
    double S_small = K_ * std::exp(x_mesh[interval_i]);
    double S_large = K_ * std::exp(x_mesh[interval_i + 1]);
    double S_larger = K_ * std::exp(x_mesh[interval_i + 2]);
    
    double V_smaller = std::exp(-a_ * (x_mesh[interval_i - 1]) - b_ * tau_final_) * u_mesh[interval_i - 1];
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_larger = std::exp(-a_ * (x_mesh[interval_i + 2]) - b_ * tau_final_) * u_mesh[interval_i + 2];
    
    // Find delta and gamma
    double delta = (V_large - V_small) / (S_large - S_small);
    double gamma = ((V_larger - V_large) / (S_larger - S_large) - (V_small - V_smaller) / (S_small - S_smaller)) / (((S_larger + S_large) / 2.) - ((S_small + S_smaller) / 2.));
    
    // THETA
    // Get dt from dtau
    double dt = 2. * dtau / (sigma_ * sigma_);
    
    // Get V at t = dt
    double V_small_prev = std::exp(-a_ * (x_mesh[interval_i]) - b_ * (tau_final_ - dtau)) * u_mesh_prev[interval_i];
    double V_large_prev = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * (tau_final_ - dtau)) * u_mesh_prev[interval_i + 1];
    double V_approx_prev = ((S_large - S0_) * V_small_prev + (S0_ - S_small) * V_large_prev) / (S_large - S_small);
    double theta = (V_approx_prev - V_approx) / dt;
    
    return std::vector<double>({delta, gamma, theta});
}

double FiniteDifferencePricer::VarianceReduction(double V_approx, std::size_t M, const Euler& euler) const {
    
    // Get corresponding European option value for BS and FD
    FiniteDifferencePricer FDPricer(S0_, K_, T_, sigma_, r_, q_);
    EuropeanOption BSPricer(0., S0_, K_, T_, sigma_, r_, q_);
    
    double FDEuro = FDPricer.Price(M, euler, European).front();
    double BSEuro = BSPricer.Put();
    
    // Adjust approximation by the pointwise difference
    return V_approx + BSEuro - FDEuro;
}

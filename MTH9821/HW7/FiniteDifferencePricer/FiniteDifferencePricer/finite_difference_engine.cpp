//
//  finite_difference_engine.cpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/12/22.
//

#include "finite_difference_engine.hpp"

#include "Decomposer.hpp"
#include "IterativeSolver.hpp"
#include "LinearSolver.hpp"

#include <iostream>

FiniteDifferenceEngine::FiniteDifferenceEngine(const FiniteDifferencePricer* pricer) : pricer_(pricer) {}

std::vector<std::vector<double>> FiniteDifferenceEngine::RunScheme(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau, const Euler& euler, const OptionType& option_type) const {
    switch (euler) {
        case Explicit:
            if (option_type == American) {
                return FullExplAmer(alpha, x_mesh, M, dtau);
            } else {
                return FullExplEuro(alpha, x_mesh, M, dtau);
            }
            break;
        case Implicit:
            if (option_type == American) {
                // !!!: Not implemented, return euro version instead
                return FullImplEuro(alpha, x_mesh, M, dtau);
            } else {
                return FullImplEuro(alpha, x_mesh, M, dtau);
            }
            break;
        case ImEx:
            if (option_type == American) {
                return FullImExAmer(alpha, x_mesh, M, dtau);
            } else {
                return FullImExEuro(alpha, x_mesh, M, dtau);
            }
            break;
    }
}

// MARK: Full run of a scheme
std::vector<std::vector<double>> FiniteDifferenceEngine::FullExplEuro(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), pricer_->boundary_tau_0_);
    
    // Print u_mesh
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->StepExplEuro(curr_tau, alpha, x_mesh, u_mesh);
        // Print u_mesh
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepExplEuro(pricer_->tau_final_, alpha, x_mesh, u_mesh);
    // Print u_mesh
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
    
}

std::vector<std::vector<double>> FiniteDifferenceEngine::FullImplEuro(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), pricer_->boundary_tau_0_);
    *(u_mesh.begin()) = pricer_->boundary_x_l_(0.);
    *(u_mesh.rbegin()) = pricer_->boundary_x_r_(0.);
    
    // Build matrix A
    mat A(mat::Zero(x_mesh.size() - 2, x_mesh.size() - 2));   // Initialize with zero matrix
    A(0, 0) = 1. + 2. * alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha;
        A(i, i - 1) = -alpha;
        A(i, i) = 1. + 2. * alpha;
    }
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->StepImplEuro(curr_tau, alpha, x_mesh, u_mesh, L, U);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepImplEuro(pricer_->tau_final_, alpha, x_mesh, u_mesh, L, U);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferenceEngine::FullImExEuro(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), pricer_->boundary_tau_0_);
    
//    this->PrintVector(u_mesh);
    
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
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->StepImExEuro(curr_tau, alpha, x_mesh, u_mesh, L, U, b_multiplier);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepImExEuro(pricer_->tau_final_, alpha, x_mesh, u_mesh, L, U, b_multiplier);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferenceEngine::FullExplAmer(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), pricer_->boundary_tau_0_);
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        // this->StepExplAmer(curr_tau, alpha, x_mesh, u_mesh);
        this->StepExplAmerFindEarlyExBoundary(curr_tau, alpha, x_mesh, u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepExplAmer(pricer_->tau_final_, alpha, x_mesh, u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

//    std::vector<std::vector<double>> FiniteDifferenceEngine::FullImplAmer(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;

std::vector<std::vector<double>> FiniteDifferenceEngine::FullImExAmer(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), pricer_->boundary_tau_0_);
//    this->PrintVector(u_mesh);
    
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
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->StepImExAmer(curr_tau, alpha, x_mesh, u_mesh, A, b_multiplier);
//        this->PrintVector(u_mesh);
        if (i % 100 == 0) {
            std::cout << i << std::endl;
        }
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->StepImExAmer(pricer_->tau_final_, alpha, x_mesh, u_mesh, A, b_multiplier);
//    this->PrintVector(u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

// MARK: Step of a full run
void FiniteDifferenceEngine::StepExplEuro(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(pricer_->boundary_x_l_(tau));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        new_u_mesh.push_back(alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1]);
    }
    
    // Right boundary
    new_u_mesh.push_back(pricer_->boundary_x_r_(tau));
    
    u_mesh = std::move(new_u_mesh);
}

void FiniteDifferenceEngine::StepImplEuro(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U) const {
    vec b(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, b.begin());
    
    // Add boundary conditions
    b(0) += pricer_->boundary_x_l_(tau) * alpha;
    b(b.size() - 1) += pricer_->boundary_x_r_(tau) * alpha;
    
    // Solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // Assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), u_mesh.begin() + 1);
    
    // Add boundary condition
    *(u_mesh.begin()) = pricer_->boundary_x_l_(tau);
    *(u_mesh.rbegin()) = pricer_->boundary_x_r_(tau);
}

void FiniteDifferenceEngine::StepImExEuro(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U, const mat& b_multiplier) const {
    // Prepare b
    vec b(u_mesh.size() - 2);
    
    // [Some matrix] * u
    vec u(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, u.begin());
    
    b = b_multiplier * u;
    
    // Boundary conditions
    b(0) += pricer_->boundary_x_l_(tau) * alpha * .5;
    b(b.size() - 1) += pricer_->boundary_x_r_(tau) * alpha * .5;
    b(0) += u_mesh.front() * alpha * .5;
    b(b.size() - 1) += u_mesh.back() * alpha * .5;
    
    // Solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // Assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), u_mesh.begin() + 1);
    
    // Add boundary condition
    *(u_mesh.begin()) = pricer_->boundary_x_l_(tau);
    *(u_mesh.rbegin()) = pricer_->boundary_x_r_(tau);
}

void FiniteDifferenceEngine::StepExplAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(pricer_->boundary_x_l_(tau));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        // Get the corresponding European option's value
        double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
        
        // Find early exercise
        double early_ex_premium = 0.;
        
        if (x_mesh[pos] < 0.) {
            early_ex_premium = pricer_->K_ * std::exp(pricer_->a_ * x_mesh[pos] + pricer_->b_ * tau) * (1. - std::exp(x_mesh[pos]));
        }
        
        // Compare and add to mesh
        
        new_u_mesh.push_back(std::max(euro_val, early_ex_premium));
    }
    
    // Right boundary
    new_u_mesh.push_back(pricer_->boundary_x_r_(tau));
    
    u_mesh = std::move(new_u_mesh);
}

void FiniteDifferenceEngine::StepExplAmerFindEarlyExBoundary(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(pricer_->boundary_x_l_(tau));
    
    std::size_t early_ex_pos = 0; // MODIFIED
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        // Get the corresponding European option's value
        double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
        
        // Find early exercise
        double early_ex_premium = 0.;
        
        if (x_mesh[pos] < 0.) {
            early_ex_premium = pricer_->K_ * std::exp(pricer_->a_ * x_mesh[pos] + pricer_->b_ * tau) * (1. - std::exp(x_mesh[pos]));
        }
        
        // Compare and add to mesh
        if (early_ex_premium > euro_val) {
            early_ex_pos = pos; // MODIFIED
        }
        
        new_u_mesh.push_back(std::max(euro_val, early_ex_premium));
    }
    
    // Right boundary
    new_u_mesh.push_back(pricer_->boundary_x_r_(tau));
    
    double S_boundary_lower = pricer_->K_ * std::exp(x_mesh[early_ex_pos]); // MODIFIED
    double S_boundary_upper = pricer_->K_ * std::exp(x_mesh[early_ex_pos + 1]); // MODIFIED
    
    std::cout << (S_boundary_lower + S_boundary_upper) / 2. << std::endl; // MODIFIED
    
    u_mesh = std::move(new_u_mesh);
}

//    void FiniteDifferenceEngine::StepImplAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U) const {}

void FiniteDifferenceEngine::StepImExAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& A, const mat& b_multiplier) const {
    // Prepare b
    vec b(u_mesh.size() - 2);
    
    // [Some matrix] * u
    vec u(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, u.begin());
    
    b = b_multiplier * u;
    
    // Boundary conditions
    b(0) += pricer_->boundary_x_l_(tau) * alpha * .5;
    b(b.size() - 1) += pricer_->boundary_x_r_(tau) * alpha * .5;
    b(0) += u_mesh.front() * alpha * .5;
    b(b.size() - 1) += u_mesh.back() * alpha * .5;
    
//    double old_right_boundary = *(u_mesh.rbegin());
    
    // Add boundary condition
    u_mesh.front() = pricer_->boundary_x_l_(tau);
    u_mesh.back() = pricer_->boundary_x_r_(tau);
    
    auto find_early_ex = [=](double x)->double {
        if (x < 0.) {
            return pricer_->K_ * std::exp(pricer_->a_ * x + pricer_->b_ * tau) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
//    vec early_ex_premium(x_mesh.size() - 2);
//    std::copy(x_mesh.cbegin() + 1, x_mesh.cend() - 1, early_ex_premium.begin());
//    std::transform(early_ex_premium.begin(), early_ex_premium.end(), early_ex_premium.begin(), find_early_ex);
//
//    // Solve linear system
//    double tolerance = std::pow(10, -6);
//    double omega = 1.2;
//
//    IterativeSolver solver(A, b, early_ex_premium);
//    vec sol(b.size());
//    std::tie(sol, std::ignore) = solver.SORProjected_lower(omega, StoppingCriterion::consecutive, tolerance, early_ex_premium);
    
    vec early_ex_premium(x_mesh.size() - 2);
    std::transform(x_mesh.cbegin() + 1, x_mesh.cend() - 1, early_ex_premium.begin(), find_early_ex);

    // Solve linear system
    double tolerance = std::pow(10, -6);
    double omega = 1.1;

    IterativeSolver solver(A, b, early_ex_premium);
    vec sol(b.size());

    std::tie(sol, std::ignore) = solver.SORProjected_lowerelementwise(omega, StoppingCriterion::consecutive, tolerance, early_ex_premium, alpha);
    
    
    // Assign to new u_mesh
    std::copy(sol.cbegin(), sol.cend(), u_mesh.begin() + 1);
}

void FiniteDifferenceEngine::PrintVector(const std::vector<double>& vec) const {
    for (double elem : vec) {
        std::cout << elem << '\t';
    }
    std::cout << std::endl;
}

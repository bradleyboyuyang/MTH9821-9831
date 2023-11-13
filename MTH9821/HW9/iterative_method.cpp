/**
 * @file linear_solver.cpp
 * @author Boyu Yang
 * @brief Jacobi, Gauss-Seidel, and SOR iterative methods implementation
 * @version 0.1
 * @date 2023-11-12
 */

#include <Eigen/Dense>
#include <tuple>
#include "iterative_method.hpp"
#include "linear_solver.hpp"

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef unsigned int uint;

// overload to use default initial guess
std::tuple<vec, uint> jacobi(const mat &A, const vec &b, const double tolerance, const StoppingCriterion criterion) {
    return jacobi(A, b, vec::Zero(A.rows()), tolerance, criterion);
}

std::tuple<vec, uint> jacobi(const mat &A, const vec &b, const vec &x_0, const double tolerance, const StoppingCriterion criterion) {
    mat D = A.diagonal().asDiagonal();
    mat R = D - A;
    mat T = D.inverse() * R;
    vec c = D.inverse() * b;
    vec x = x_0, x_old;
    uint iter = 0;
    // loop and stop based on the stopping criterion
    if (criterion == StoppingCriterion::consecutive) {
        do {
            x_old = x;
            x = T * x + c;
            ++iter;
        } while ((x - x_old).norm() > tolerance);
    } else if (criterion == StoppingCriterion::residual) {
        do {
            x_old = x;
            x = T * x + c;
            ++iter;
        } while ((A * x - b).norm() > tolerance);
    }

    return std::make_tuple(x, iter);
}

std::tuple<vec, uint> gauss_seidel(const mat &A, const vec &b, const double tolerance, const StoppingCriterion criterion) {
    return gauss_seidel(A, b, vec::Zero(A.rows()), tolerance, criterion);
}

std::tuple<vec, uint> gauss_seidel(const mat &A, const vec &b, const vec &x_0, const double tolerance, const StoppingCriterion criterion) {
    // A = D + L + U
    mat D = A.diagonal().asDiagonal();
    mat L = mat::Zero(A.rows(), A.cols());
    mat U = mat::Zero(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            if (j < i) {
                L(i, j) = A(i, j);
            } else if (j > i) {
                U(i, j) = A(i, j);
            }
        }
    }
    // calculate inverse or use forward substitution
    // mat T = -(D + L).inverse() * U;
    // vec c = (D + L).inverse() * b;
    mat T;
    vec c = forward_subst(D + L, b);

    vec x = x_0, x_old;
    uint iter = 0;
    // loop and stop based on the stopping criterion
    if (criterion == StoppingCriterion::consecutive) {
        do {
            x_old = x;
            x = -forward_subst(D + L, U*x) + c;
            // x = T * x + c;
            ++iter;
        } while ((x - x_old).norm() > tolerance);
    } else if (criterion == StoppingCriterion::residual) {
        do {
            x_old = x;
            x = -forward_subst(D + L, U*x) + c;
            // x = T * x + c;
            ++iter;
        } while ((A * x - b).norm() > tolerance);
    }

    return std::make_tuple(x, iter);
}


std::tuple<vec, uint> sor(const mat &A, const vec &b, const double omega, const double tolerance, const StoppingCriterion criterion) {
    return sor(A, b, vec::Zero(A.rows()), omega, tolerance, criterion);
}

std::tuple<vec, uint> sor(const mat &A, const vec &b, const vec &x_0, const double omega, const double tolerance, const StoppingCriterion criterion) {
    // A = D + L + U
    mat D = A.diagonal().asDiagonal();
    mat L = mat::Zero(A.rows(), A.cols());
    mat U = mat::Zero(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            if (j < i) {
                L(i, j) = A(i, j);
            } else if (j > i) {
                U(i, j) = A(i, j);
            }
        }
    }
    // calculate inverse or use forward substitution
    // mat T = -(D + omega*L).inverse() * ((1-omega)*D - omega*U);
    // vec c = omega * (D + omega*L).inverse() * b;
    mat T;
    vec c = omega * forward_subst(D + omega*L, b);

    vec x = x_0, x_old;
    uint iter = 0;
    // loop and stop based on the stopping criterion
    if (criterion == StoppingCriterion::consecutive) {
        do {
            x_old = x;
            x = forward_subst(D + omega*L, ((1-omega)*D-omega*U)*x) + c;
            // x = T * x + c;
            ++iter;
        } while ((x - x_old).norm() > tolerance);
    } else if (criterion == StoppingCriterion::residual) {
        do {
            x_old = x;
            x = forward_subst(D + omega*L, ((1-omega)*D-omega*U)*x) + c;
            // x = T * x + c;
            ++iter;
        } while ((A * x - b).norm() > tolerance);
    }

    return std::make_tuple(x, iter);
}



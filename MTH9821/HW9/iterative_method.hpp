/**
 * @file linear_solver.cpp
 * @author Boyu Yang
 * @brief Header file for iterative methods
 * @version 0.1
 * @date 2023-11-12
 */

#ifndef iterative_methods_hpp
#define iterative_methods_hpp

#include <Eigen/Dense>
#include <tuple>
#include <iostream>

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef unsigned int uint;
enum class StoppingCriterion { consecutive, residual };

// overload to use default initial guess
std::tuple<vec, uint> jacobi(const mat &A, const vec &b, const double tolerance, const StoppingCriterion criterion);
std::tuple<vec, uint> jacobi(const mat &A, const vec &b, const vec &x_0, const double tolerance, const StoppingCriterion criterion);
std::tuple<vec, uint> gauss_seidel(const mat &A, const vec &b, const double tolerance, const StoppingCriterion criterion);
std::tuple<vec, uint> gauss_seidel(const mat &A, const vec &b, const vec &x_0, const double tolerance, const StoppingCriterion criterion);
std::tuple<vec, uint> sor(const mat &A, const vec &b, const double omega, const double tolerance, const StoppingCriterion criterion);
std::tuple<vec, uint> sor(const mat &A, const vec &b, const vec &x_0, const double omega, const double tolerance, const StoppingCriterion criterion);

#endif /* iterative_methods_hpp */

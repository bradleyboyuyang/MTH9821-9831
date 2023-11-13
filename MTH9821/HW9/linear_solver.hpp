/**
 * @file linear_solver.hpp
 * @author Yicheng Sun
 * @brief header file for linear solver
 * @version 0.1
 * @date 2023-11-09
 */

#ifndef linear_solver_hpp
#define linear_solver_hpp

#include <stdio.h>
#include <Eigen/Dense>

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1, -1, uint> permutation;

vec forward_subst(const mat &L, const vec &b);
vec backward_subst(const mat &U, const vec &b);

#endif /* linear_solver_hpp */

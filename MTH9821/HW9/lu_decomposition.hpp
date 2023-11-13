/**
 * @file lu_decomposition.hpp
 * @author Yicheng Sun
 * @brief header file of LU decomposition
 * @version 0.1
 * @date 2023-11-09
 */

#ifndef lu_decomposition_hpp
#define lu_decomposition_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <tuple>

typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1, -1, uint> permutation;

std::tuple<mat, mat> lu_no_pivoting(mat A);
std::tuple<permutation, mat, mat> lu_row_pivoting(mat A);

#endif /* lu_decomposition_hpp */

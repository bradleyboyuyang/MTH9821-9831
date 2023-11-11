//
//  lu_decomposition.hpp
//  nla
//
//  Created by Yicheng.Sun on 9/11/2023.
//

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

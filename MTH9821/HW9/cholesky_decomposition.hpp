/**
 * @file cholesky_decomposition.hpp
 * @author Yicheng Sun
 * @brief header file of cholesky decomposition
 * @version 0.1
 * @date 2023-11-09
 */

#ifndef cholesky_decomposition_hpp
#define cholesky_decomposition_hpp

#include <stdio.h>
#include <Eigen/Dense>

typedef Eigen::MatrixXd mat;

mat cholesky(mat A);

#endif /* cholesky_decomposition_hpp */

//
//  cholesky_decomposition.hpp
//  nla
//
//  Created by Yicheng.Sun on 9/11/2023.
//

#ifndef cholesky_decomposition_hpp
#define cholesky_decomposition_hpp

#include <stdio.h>
#include <Eigen/Dense>

typedef Eigen::MatrixXd mat;

mat cholesky(mat A);

#endif /* cholesky_decomposition_hpp */

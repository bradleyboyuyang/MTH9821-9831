//
//  Decomposer.hpp
//  MTH9821_HW9_Group6
//
//  Created by 王明森 on 11/10/22.
//

#ifndef Decomposer_hpp
#define Decomposer_hpp

#include "EigenCommonHeader.h"
#include <tuple>

class Decomposer {
public:
    std::tuple<mat, mat> lu_no_pivoting (const mat& A) const;
    std::tuple<mat, mat, mat> lu_row_pivoting (const mat& A) const;
    mat cholesky(const mat& A) const;
};

#endif /* Decomposer_hpp */

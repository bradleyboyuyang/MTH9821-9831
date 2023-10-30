//
//  Decomposer.cpp
//  MTH9821_HW9_Group6
//
//  Created by 王明森 on 11/10/22.
//

#include "Decomposer.hpp"

std::tuple<mat, mat> Decomposer::lu_no_pivoting (const mat& A) const {
    
    mat A_copy = A;
    
    mat L = A.triangularView<Eigen::UnitLower>();
    mat U = A.triangularView<Eigen::Upper>();
    
    for (std::size_t i = 0; i < A.rows() - 1; i++) {
        for (std::size_t k = i; k < A.rows(); k++) {
            U(i, k) = A_copy(i, k);
            L(k, i) = A_copy(k, i) / U(i, i);
        }
        for (std::size_t r = i + 1; r < A.rows(); r++) {
            for (std::size_t c = i + 1; c < A.rows(); c++) {
                A_copy(r, c) -= L(r, i) * U(i, c);
            }
        }
        U(A.rows() - 1, A.rows() - 1) = A_copy(A.rows() - 1, A.rows() - 1);
    }
    
    return std::make_tuple(L, U);
}

std::tuple<mat, mat, mat> Decomposer::lu_row_pivoting (const mat& A) const {
    Eigen::PartialPivLU<mat> decomposer(A);
    
    // Get P and LU
    mat P = decomposer.permutationP();
    mat LU = decomposer.matrixLU();
    
    // Recover L and U
    mat L = LU.triangularView<Eigen::UnitLower>();
    mat U = LU.triangularView<Eigen::Upper>();
    
    return std::make_tuple(P, L, U);
}

mat Decomposer::cholesky(const mat& A) const {
    // A = UtU
    return A.llt().matrixL().transpose();
}

//
//  lu_decomposition.cpp
//  nla
//
//  Created by Yicheng.Sun on 9/11/2023.
//

#include "lu_decomposition.hpp"
#include <iostream>

std::tuple<mat, mat> lu_no_pivoting(mat A) {
    int n = static_cast<int>(A.rows());
    mat L = mat::Identity(n, n);
    mat U = A;

    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            L(i, k) = U(i, k) / U(k, k);
            for (int j = k; j < n; ++j) {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }
    }

    return std::make_tuple(L, U);
}


std::tuple<permutation, mat, mat> lu_row_pivoting(mat A) {
    int n = static_cast<int>(A.rows());
    
    permutation P(n);
    P.setIdentity();
    mat L = mat::Identity(n, n);
    mat U = A;

    for (int k = 0; k < n - 1; ++k) {
        // Partial Pivoting
        int pivot_index = k;
        double max_val = std::abs(A(k, k));

        for (int j = k + 1; j < n; j++) {
            if (std::abs(U(j, k)) > max_val) {
                max_val = std::abs(U(j, k));
                pivot_index = j;
            }
        }

        if (pivot_index != k) {
            // Swap rows in U
            U.row(k).swap(U.row(pivot_index));

            // Swap rows in L
            L.block(k, 0, 1, k).swap(L.block(pivot_index, 0, 1, k));
            
            int tmp = P.indices()[k];
            P.indices()[k] = P.indices()[pivot_index];
            P.indices()[pivot_index] = tmp;
        }
    
        // Gaussian elimination
        for (int i = k + 1; i < n; ++i) {
            L(i, k) = U(i, k) / U(k, k);
            U.row(i) -= L(i, k) * U.row(k);
        }
    }

    return std::make_tuple(P, L, U);
}

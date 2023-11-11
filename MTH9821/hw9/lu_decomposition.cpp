//
//  lu_decomposition.cpp
//  nla
//
//  Created by Yicheng.Sun on 9/11/2023.
//

#include "lu_decomposition.hpp"

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


std::tuple<permutation, mat, mat> lu_row_pivoting(const mat& A) {
    int n = static_cast<int>(A.rows());
    permutation P(n);
    mat L = mat::Identity(n, n);
    mat U = A;

    for (int k = 0; k < n - 1; ++k) {
        // Partial Pivoting
        int pivot_index = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(U(i, k)) > std::abs(U(pivot_index, k))) {
                pivot_index = i;
            }
        }

        if (pivot_index != k) {
            // Swap rows in U
            U.row(k).swap(U.row(pivot_index));

            // Swap rows in L
            L.row(k).swap(L.row(pivot_index));

            // Apply permutation to P
            P.applyTranspositionOnTheLeft(k, pivot_index);
        }

        // Gaussian elimination
        for (int i = k + 1; i < n; ++i) {
            L(i, k) = U(i, k) / U(k, k);
            U.row(i) -= L(i, k) * U.row(k);
        }
    }

    return std::make_tuple(P, L, U);
}

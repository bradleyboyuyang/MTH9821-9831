//
//  cholesky_decomposition.cpp
//  nla
//
//  Created by Yicheng.Sun on 9/11/2023.
//

#include "cholesky_decomposition.hpp"

mat cholesky(mat A) {
    int n = static_cast<int>(A.rows());
    mat L(n, n);

    for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < j; ++k) {
            sum += L(j, k) * L(j, k);
        }

        L(j, j) = std::sqrt(A(j, j) - sum);

        for (int i = j + 1; i < n; ++i) {
            sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += L(i, k) * L(j, k);
            }
            L(i, j) = (A(i, j) - sum) / L(j, j);
        }
    }

    return L;
}


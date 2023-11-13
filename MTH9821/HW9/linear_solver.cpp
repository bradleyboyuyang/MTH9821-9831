/**
 * @file linear_solver.cpp
 * @author Yicheng Sun
 * @brief implementation file for linear solver
 * @version 0.1
 * @date 2023-11-09
 */

#include "linear_solver.hpp"

vec forward_subst(const mat &L, const vec &b) {
    int n = static_cast<int>(L.rows());
    vec y(n);

    for (int i = 0; i < n; ++i) {
        y(i) = b(i);
        for (int j = 0; j < i; ++j) {
            y(i) -= L(i, j) * y(j);
        }
        y(i) /= L(i, i);
    }

    return y;
}

vec backward_subst(const mat &U, const vec &b) {
    int n = static_cast<int>(U.rows());
    vec x(n);

    for (int i = n - 1; i >= 0; --i) {
        x(i) = b(i);
        for (int j = i + 1; j < n; ++j) {
            x(i) -= U(i, j) * x(j);
        }
        x(i) /= U(i, i);
    }
    return x;
}

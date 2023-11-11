//
//  main.cpp
//  nla
//
//  Created by Yicheng.Sun on 9/11/2023.
//

#include <iostream>
#include <iostream>
#include "linear_solver.hpp"
#include "lu_decomposition.hpp"
#include "cholesky_decomposition.hpp"

int main() {
    // test linear_solver
//    mat L(3, 3); // Lower triangular matrix
//    mat U(3, 3); // Upper triangular matrix
//    vec b(3);    // Right-hand side vector
//
//    // Fill in L, U, and b with appropriate values
//    L << 1, 0, 0,
//         2, 3, 0,
//         4, 5, 6;
//
//    U << 1, 2, 3,
//         0, 4, 5,
//         0, 0, 6;
//
//    b << 6, 9, 15;
//
//    // Example forward substitution
//    vec y = forward_subst(L, b);
//
//    // Example backward substitution
//    vec x = backward_subst(U, b);

    // Display the results
//    std::cout << "Lower Triangular Matrix L:\n" << L << "\n";
//    std::cout << "Upper Triangular Matrix U:\n" << U << "\n";
//    std::cout << "Right-hand Side Vector b: " << b.transpose() << "\n";
//    std::cout << "Forward Substitution Result (y): " << y.transpose() << "\n";
//    std::cout << "Backward Substitution Result (x): " << x.transpose() << "\n";

    
//     test lu_decomp
    mat A(3, 3); // Input matrix
    A << 2, 1, -1,
         -3, -1, 2,
         -2, 1, 2;

    // LU decomposition without pivoting
//    auto [L_no_pivot, U_no_pivot] = lu_no_pivoting(A);

    // LU decomposition with row pivoting
    permutation P;
    mat L, U;

    std::tie(P, L, U) = lu_row_pivoting(A);

    // Output the results
    std::cout << "Permutation matrix P:\n" << P.toDenseMatrix() << "\n\n";
    std::cout << "Lower triangular matrix L:\n" << L << "\n\n";
    std::cout << "Upper triangular matrix U:\n" << U << "\n";
//    auto [P_row_pivot, L_row_pivot, U_row_pivot] = lu_row_pivoting(A);

    // Display the results
//    std::cout << "Input Matrix A:\n" << A << "\n";
//    std::cout << "LU Decomposition without Pivoting (A = LU):\n";
//    std::cout << "L:\n" << L_no_pivot << "\nU:\n" << U_no_pivot << "\n";

//    std::cout << "LU Decomposition with Row Pivoting (PA = LU):\n";
//    std::cout << "P:\n" << P_row_pivot.indices() << "\n";
//    std::cout << "L:\n" << L_row_pivot << "\nU:\n" << U_row_pivot << "\n";

    
    // test cholesky_decomp
//    mat A_new(3, 3); // SPD matrix
//    A_new << 4, 12, -16,
//        12, 37, -43,
//        -16, -43, 98;
//
//    // Cholesky decomposition
//    mat L_new = cholesky(A_new);
//
//    // Display the result
//    std::cout << "Original SPD Matrix A:\n" << A_new << "\n";
//    std::cout << "Cholesky Decomposition (L, where A = LL^T):\n" << L_new << "\n";

    
    return 0;
}

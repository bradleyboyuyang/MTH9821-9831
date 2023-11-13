/**
 * @file main.cpp
 * @author Boyu Yang, Yicheng Sun
 * @brief main file to answer all questions
 * @version 0.1
 * @date 2023-11-12
 */

#include <iostream>
#include "linear_solver.hpp"
#include "lu_decomposition.hpp"
#include "cholesky_decomposition.hpp"
#include "iterative_method.hpp"

int main() {
    /// Problem 1: test forward_subst and backward_subst
    mat L(3, 3); // Lower triangular matrix
    mat U(3, 3); // Upper triangular matrix
    vec b(3);    // Right-hand side vector

    // test case
    L << 1, 0, 0,
         2, 3, 0,
         4, 5, 6;

    U << 1, 2, 3,
         0, 4, 5,
         0, 0, 6;

    b << 6, 9, 15;

    // Example forward substitution
    vec y = forward_subst(L, b);

    // Example backward substitution
    vec x = backward_subst(U, b);

    std::cout << "===========================Problem 1===========================\n";
    std::cout << "Lower Triangular Matrix L:\n" << L << "\n";
    std::cout << "Upper Triangular Matrix U:\n" << U << "\n";
    std::cout << "Right-hand Side Vector b: " << b.transpose() << "\n";
    std::cout << "Forward Substitution Result (y): " << y.transpose() << "\n";
    std::cout << "Backward Substitution Result (x): " << x.transpose() << "\n";

    /// Problem 2: test lu_no_pivoting and lu_row_pivoting    
    mat A(3, 3); // Input matrix
    A << 4, 3, 0,
         3, 4, -1,
         0, -1, 4;
    b << 24, 30, -24;

    // LU decomposition without pivoting
    auto [L_no_pivot, U_no_pivot] = lu_no_pivoting(A);

    std::cout << "===========================Problem 2===========================\n";
    std::cout << "Input Matrix A:\n" << A << "\n";
    std::cout << "LU Decomposition without Pivoting (A = LU):\n";
    std::cout << "L:\n" << L_no_pivot << "\nU:\n" << U_no_pivot << "\n";
    // solve Ax = b
    x = backward_subst(U_no_pivot, forward_subst(L_no_pivot, b));
    std::cout << "Solution x (without pivoting): " << x.transpose() << "\n";

    auto [P_row_pivot, L_row_pivot, U_row_pivot] = lu_row_pivoting(A);
    std::cout << "LU Decomposition with Row Pivoting (PA = LU):\n";
    std::cout << "P:" << P_row_pivot.indices().transpose() << "\n";
    std::cout << "L:\n" << L_row_pivot << "\nU:\n" << U_row_pivot << "\n";
    // solve Ax = b
    x = backward_subst(U_row_pivot, forward_subst(L_row_pivot, P_row_pivot * b));
    std::cout << "Solution x (with pivoting): " << x.transpose() << "\n";

    /// Problem 3: test cholesky
    mat L3 = cholesky(A);
    std::cout << "===========================Problem 3===========================\n";
    std::cout << "Original SPD Matrix A:\n" << A << "\n";
    std::cout << "Cholesky Decomposition (L, where A = LL^T):\n" << L3 << "\n";
    // solve Ax = b
    x = backward_subst(L3.transpose(), forward_subst(L3, b));
    std::cout << "Solution x (with Cholesky): " << x.transpose() << "\n";
    // we can see the solutions are the same with LU decomposition with/without pivoting

    /// Problem 4: test jacobi, gauss_seidel, and SOR method
    // Initial guess
    vec x_0 = vec::Ones(3);
    // tolerance
    double tolerance = 1e-6;

    // Testing Jacobi Method
    auto [x_jacobi, iterations_jacobi] = jacobi(A, b, x_0, tolerance, StoppingCriterion::consecutive);
    std::cout << "===========================Problem 4===========================\n";
    std::cout << "Jacobi Method:" << std::endl;
    std::cout << "Solution: " << x_jacobi.transpose() << std::endl;
    std::cout << "Iterations: " << iterations_jacobi << std::endl;

    // Testing Gauss-Siedel Method
    auto [x_gauss, iterations_gauss] = gauss_seidel(A, b, x_0, tolerance, StoppingCriterion::consecutive);
    std::cout << "\nGauss-Siedel Method:" << std::endl;
    std::cout << "Solution: " << x_gauss.transpose() << std::endl;
    std::cout << "Iterations: " << iterations_gauss << std::endl;

    // Testing SOR Method with a relaxation factor (omega)
    double omega = 1.25;  // Example value, you can adjust it
    auto [x_sor, iterations_sor] = sor(A, b, x_0, omega, tolerance, StoppingCriterion::consecutive);
    std::cout << "\nSOR Method:" << std::endl;
    std::cout << "Solution: " << x_sor.transpose() << std::endl;
    std::cout << "Iterations: " << iterations_sor << std::endl;
    
    /// Problem 5: see the proof in the attached picture

    /// Problem 6: check performance of the three methods on a 14x14 tridiagonal matrix
    // Define the 14x14 matrix A
    mat A_new = mat::Zero(14, 14);
    for (int i = 0; i < 14; ++i) {
        A_new(i, i) = 2;
        if (i > 0) {
            A_new(i, i - 1) = -1;
            A_new(i - 1, i) = -1;
        }
    }

    // Define the vector b
    vec b_new(14);
    for (int i = 0; i < 14; ++i) {
        b_new(i) = i * i;
    }

    // Initial guess and tolerance
    vec x_0_new = vec::Constant(14, 1);
    tolerance = 1e-6;

    // Jacobi Iteration: residual-based vs consecutive approximation stopping criterion
    std::cout << "===========================Problem 6===========================\n";
    auto [x_jacobi_residual, iter_jacobi_residual] = jacobi(A_new, b_new, x_0_new, tolerance, StoppingCriterion::residual);
    std::cout << "Jacobi method (Residual): \nIterations " << iter_jacobi_residual << ", Solution " << x_jacobi_residual.transpose() << "\n";
    auto [x_jacobi_consecutive, iter_jacobi_consecutive] = jacobi(A_new, b_new, x_0_new, tolerance, StoppingCriterion::consecutive); 
    std::cout << "Jacobi method (Consecutive): \nIterations " << iter_jacobi_consecutive << ", Solution " << x_jacobi_consecutive.transpose() << "\n";

    // Gauss-Siedel Iteration: residual-based vs consecutive approximation stopping criterion
    auto [x_gauss_residual, iter_gauss_residual] = gauss_seidel(A_new, b_new, x_0_new, tolerance, StoppingCriterion::residual); 
    std::cout << "Gauss-Siedel (Residual): \nIterations " << iter_gauss_residual << ", Solution " << x_gauss_residual.transpose() << "\n";

    auto [x_gauss_consecutive, iter_gauss_consecutive] = gauss_seidel(A_new, b_new, x_0_new, tolerance, StoppingCriterion::consecutive); 
    std::cout << "Gauss-Siedel (Consecutive): \nIterations " << iter_gauss_consecutive << ", Solution " << x_gauss_consecutive.transpose() << "\n";

    // SOR Iteration: residual-based vs consecutive approximation stopping criterion
    omega = 1.15;
    auto [x_sor_residual, iter_sor_residual] = sor(A_new, b_new, x_0_new, omega, tolerance, StoppingCriterion::residual);
    std::cout << "SOR (Residual): \nIterations " << iter_sor_residual << ", Solution " << x_sor_residual.transpose() << "\n";
    auto [x_sor_consecutive, iter_sor_consecutive] = sor(A_new, b_new, x_0_new, omega, tolerance, StoppingCriterion::consecutive);
    std::cout << "SOR (Consecutive): \nIterations " << iter_sor_consecutive << ", Solution " << x_sor_consecutive.transpose() << "\n";

    // residual based, solve Ax=b using SOR for omega = 1.02 : .02 : 1.98, report iterations to convergence
    std::cout << "SOR iteration with w = 1.02 : .02 : 1.98:" << std::endl;
    for (double omega = 1.02; omega < 2; omega += 0.02) {
        auto [x_sor_residual, iter_sor_residual] = sor(A_new, b_new, x_0_new, omega, tolerance, StoppingCriterion::residual);
        std::cout << omega << ", " << iter_sor_residual << std::endl;
    }
    // Comment: w is a relaxation factor, which can be used to reduce the norm of residual vector and accelerate the convergence of the SOR method.
    // When conducting over-relaxation (w > 1), the SOR method converges faster as w increases. However, when w is too large, 
    // The solution x may face overshooting problem, making the iteration change too fast and again require more iterations to converge.
    // Overall, the optimal w is around 1.66, which only requires 54 iterations to converge.
    return 0;
}

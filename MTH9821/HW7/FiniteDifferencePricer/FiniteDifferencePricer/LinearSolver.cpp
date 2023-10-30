//
//  LinearSolver.cpp
//  MTH9821_HW9_Group6
//
//  Created by 王明森 on 11/10/22.
//

#include "LinearSolver.hpp"
#include <exception>

vec LinearSolver::ForwardSub(const mat& L, const vec& b) const {
    assert (L.isLowerTriangular());
    return L.triangularView<Eigen::Lower>().solve(b);
}

vec LinearSolver::BackwardSub(const mat& U, const vec& b) const {
    assert (U.isUpperTriangular());
    return U.triangularView<Eigen::Upper>().solve(b);
}

vec LinearSolver::LUSolve(const mat& A, const vec& b) const {
    return A.partialPivLu().solve(b);
}

vec LinearSolver::CholeskySolve(const mat& A, const vec& b) const {
    return A.llt().solve(b);
}

vec LinearSolver::DiagonalSolve_vec(const mat& D, const vec& b) const {
//    assert (D.isDiagonal());
    
    // Construct arrays for coefficient-wise operation
    Eigen::ArrayXd b_arr(b);
    Eigen::ArrayXd D_arr(D.diagonal());
    
    return vec(b_arr / D_arr);
}

mat LinearSolver::DiagonalSolve(const mat& D, const mat& B) const {
//     assert (D.isDiagonal());
    
    mat res(B);
    
    for (std::size_t i = 0; i < B.cols(); i++) {
        res(Eigen::all, i) = this->DiagonalSolve_vec(D, B(Eigen::all, i));
    }
    
    return res;
}

mat LinearSolver::ForwardSub_mat(const mat& L, const mat& B) const {
    mat res(B);
    
    for (std::size_t i = 0; i < B.cols(); i++) {
        res(Eigen::all, i) = this->ForwardSub(L, B(Eigen::all, i));
    }
    
    return res;
}

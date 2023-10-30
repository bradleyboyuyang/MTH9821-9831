//
//  LinearSolver.hpp
//  MTH9821_HW9_Group6
//
//  Created by 王明森 on 11/10/22.
//

#ifndef LinearSolver_hpp
#define LinearSolver_hpp

#include "EigenCommonHeader.h"

class LinearSolver {
public:
    vec ForwardSub(const mat& L, const vec& b) const;
    vec BackwardSub(const mat& U, const vec& b) const;
    vec LUSolve(const mat& A, const vec& b) const;
    vec CholeskySolve(const mat& A, const vec& b) const;
    
    vec DiagonalSolve_vec(const mat& D, const vec& b) const;
    mat DiagonalSolve(const mat& D, const mat& B) const;
    
    mat ForwardSub_mat(const mat& L, const mat& B) const;
};

#endif /* LinearSolver_hpp */

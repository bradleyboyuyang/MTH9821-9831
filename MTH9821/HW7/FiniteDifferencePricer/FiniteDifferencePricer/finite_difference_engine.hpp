//
//  finite_difference_engine.hpp
//  FiniteDifferencePricer
//
//  Created by 王明森 on 12/12/22.
//

#ifndef finite_difference_engine_hpp
#define finite_difference_engine_hpp

#include <stdio.h>
#include <vector>
#include "finite_difference_pricer.hpp"
#include "EigenCommonHeader.h"

class FiniteDifferenceEngine {
private:
    const FiniteDifferencePricer* pricer_;
    
public:
    FiniteDifferenceEngine(const FiniteDifferencePricer* pricer);
    ~FiniteDifferenceEngine() = default;
    
    std::vector<std::vector<double>> RunScheme(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau, const Euler& euler, const OptionType& option_type) const;
    
    void PrintVector(const std::vector<double>& vec) const;
    
private:
    // Full run of a scheme
    std::vector<std::vector<double>> FullExplEuro(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FullImplEuro(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FullImExEuro(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FullExplAmer(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
//    std::vector<std::vector<double>> FullImplAmer(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FullImExAmer(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    // Step of a full run
    void StepExplEuro(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    
    void StepImplEuro(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U) const;
    
    void StepImExEuro(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U, const mat& b_multiplier) const;
    
    void StepExplAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    
//    void StepImplAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U) const;
    
    void StepImExAmer(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& A, const mat& b_multiplier) const;
    
    void StepExplAmerFindEarlyExBoundary(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    
};

#endif /* finite_difference_engine_hpp */

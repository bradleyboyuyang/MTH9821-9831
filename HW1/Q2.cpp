#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>

#include <iostream>
#include <cmath>
#include <limits>



// Constants
const double S0 = 39.0;
const double K = 39.0;
const double B = 36;
const double T = 0.75;
const double sigma = 0.25;
const double q = 0.01;
const double r = 0.02;

#include <iostream>
#include <cmath>
#include <limits>

// Constants for approximations to inverse normal from Figure 2.12
const double a0 = 2.50662823884;
const double a1 = -18.61500062529;
const double a2 = 41.39119773534;
const double a3 = -25.44106049637;
const double b0 = -8.47351093090;
const double b1 = 23.08336743743;
const double b2 = -21.06224101826;
const double b3 = 3.13082909833;
const double c0 = 0.3374754822726147;
const double c1 = 0.9761690190917186;
const double c2 = 0.1607979714918209;
const double c3 = 0.0276438810333863;
const double c4 = 0.0038405729373609;
const double c5 = 0.0003951896511919;
const double c6 = 0.0000321767881768;
const double c7 = 0.0000002888167364;
const double c8 = 0.0000003960315187;

// Function to compute the inverse of the cumulative distribution function of a standard normal distribution
double norm_ppf(double u) {
    double y = u - 0.5;
    double r, x;

    if (std::abs(y) < 0.42) {
        r = y * y;
        x = y * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.0);
    } else {
        r = u;
        if (y > 0) r = 1 - u;
        r = std::log(-std::log(r));
        x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
        if (y < 0) x = -x;
    }

    
    return x;
}


// Function to calculate the cumulative distribution function of the standard normal distribution
double norm_cdf(double x) {
    return 0.5 * erfc(-x / sqrt(2.0));
}
// Function to calculate the Black-Scholes formula for European call option
double C(double S, double K) {
    double d1 = (log(S / K) + (r - q + sigma * sigma / 2) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    return S * exp(-q * T) * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
}

// Function to calculate the option value V
double V(double S, double K) {
    double a = (r - q) / (sigma * sigma) - 0.5;
    return C(S, K) - pow((B / S), (2 * a)) * C((B * B) / S, K);
}

// Function to calculate the option value using Monte Carlo simulation
double V_mc(const std::vector<double>& path, double K, double B) {
    
    if (*std::min_element(path.begin(), path.end()) > B) {
        double V_mc = std::max(0.0, path.back() - K) * exp(-r * T);
        // std::cout << V_mc;
        return V_mc;
    }
    else {
        return 0.0;
    }
}

// Function to calculate the option values for down-and-out calls
std::vector<double> valuing_down_and_out_calls(int n, int m) {
    int N = n * m;
    double dt = T / m;
    std::vector<std::vector<double>> S(n, std::vector<double>(m, S0));

    for (int i = 0; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            double Z = norm_ppf(rand() / (RAND_MAX + 1.0));
            S[i][j] = S[i][j - 1] * exp((r - q - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * Z);
        }
    }
    double V_calc = V(S[0][0], K);
    double V_est = 0.0;
    for (int i = 0; i < n; ++i) {
        
        V_est += V_mc(S[i], K, B);
    }
    V_est /= n;

    double approximation_error = fabs(V_est - V_calc);
    return {V_calc, V_est, approximation_error};
}

int main() {
    for (int k = 0; k <= 9; ++k) {
        int n = 50 * pow(2, k);
        int m = 200;
        std::vector<double> results = valuing_down_and_out_calls(n, m);
        std::cout << n * m << " " << m << " " << n << " " << results[0] << " " << results[1] << " " << results[2] << std::endl;
    }

    for (int k = 0; k <= 9; ++k) {
        int Nk = 10000 * pow(2, k);
        int mk = ceil(pow(Nk, 1.0 / 3) * pow(T, 2.0 / 3));
        int nk = floor(Nk / mk);
        std::vector<double> results = valuing_down_and_out_calls(nk, mk);
        std::cout << Nk << " " << mk << " " << nk << " " << results[0] << " " << results[1] << " " << results[2] << std::endl;
    }

    return 0;
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include "payoff.h"
#include "option_value.h"
#include "trinomial.h"

// print a vector of doubles
void print(const std::vector<double>& v)
{
	for (int i = 0; i < v.size(); i++)
	{
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}

int main()
{
	// parameters
	double T = 7.0 / 12.0;
	double K = 40;
	double S = 42;
	double r = 0.04;
	double q = 0.015;
	double vol = 0.28;
	double B = 36;
	// get theoretical value
	double a = (r-q)/(vol*vol)-0.5;
	OptionValue exact = std::get<0>(BlackScholes(T, K, S, r, q, vol,0.0));
	OptionValue modified = std::get<0>(BlackScholes(T, K, B*B/S, r, q, vol, 0.0));
	double exactval = exact.price - std::pow(B/S, 2*a)*modified.price;
	std::cout << "exact value: " << exactval << std::endl;

	// print precision
	int p = 6;

	CallPayoff barrierCall(K);
	TrinomialTree t1(barrierCall, T, S, r, q, vol);

	std::cout << std::fixed
		<< std::setprecision(p);

	std::ofstream ofs("output.csv"); // Open the output file

	//find optimal N for binomial tree
	std::vector<double> K_value;
	for (int i = 10; i <= 1000; i++)
    {
        double value = t1.barrierOption_trinomial_tree(i, B).price;
        ofs << i << ","  << value << "," << std::abs(value-exactval) << std::endl;
    }

	ofs.close(); // Close the output file

	// calculate optimal N for k = 2-13
	for (int k = 2; k <= 13; k++){
		int N = floor(k*k*(3*vol*vol*T)/(std::log(S/B)*std::log(S/B)));
		std::cout << "K=" << k << ", N=" << N << ", " << "option value:" << t1.barrierOption_trinomial_tree(N, B).price << std::endl;
	}


	return 0;
}

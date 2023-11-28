#pragma once
#include<vector>

enum Type { European, American };
enum CP { call, put };

struct Option
{
	Type type;
	CP cp;
	double S0;
	double K;
	double T;
	double sigma;
	double r;
	std::vector<std::vector<double>> div_vec;
};

std::vector<double> european_put_fixed_depth_binomial_tree(Option& op, int N);
std::vector<double> european_put_binomial_tree(Option& op, int N0, double tol);

std::vector<double> european_call_fixed_depth_binomial_tree(Option& op, int N);
std::vector<double> european_call_binomial_tree(Option& op, int N0, double tol);

std::vector<double> american_put_fixed_depth_binomial_tree(Option& op, int N);
std::vector<double> american_put_binomial_tree(Option& op, int N0, double tol);

std::vector<double> american_call_fixed_depth_binomial_tree(Option& op, int N);
std::vector<double> american_call_binomial_tree(Option& op, int N0, double tol);

std::vector<double> BinomialTree(Option& op, int N0, double tol);
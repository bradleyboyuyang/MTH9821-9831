#include<cmath>
#include<algorithm>
#include "binomial_tree.h"

// fixed depth binomial tree pricer for european call option with discrete dividend payment
std::vector<double> european_call_fixed_depth_binomial_tree(Option& op, int N)
{
	double dt = op.T / N;
	double disc = std::exp(-op.r * dt);
	double u = std::exp(op.sigma * std::sqrt(dt));
	double d = 1 / u;
	double pu = (std::exp(op.r*dt) - d) / (u - d);
	double pd = 1 - pu;

	int num_div = op.div_vec.size(); // the total number of dividend payments

	// store the first six option prices
	std::vector<std::vector<double>> option_result_vec = { {0},{0,0},{0,0,0} };
	std::vector<double> V(N + 1);

	// store the first six stock prices
	std::vector<std::vector<double>> stock_result_vec = { {0},{0,0},{0,0,0} };
	// An intermediate vector used for computing stock_result_vec
	std::vector<double> s_div_vec(6);

	// Calculate the corrected initial stock price without dividend
	double s_no_div0 = op.S0;
	for (int i = 0; i < num_div; i++)
	{
		if (op.div_vec[i][0] == 0)
		{
			s_no_div0 *= (1 - op.div_vec[i][1]);
		}
		else
		{
			s_no_div0 -= op.div_vec[i][1] * std::exp(-op.r * op.div_vec[i][2]);
		}
	}

	// Initialize V with option payment at maturity T coresponding the non-dividend-paying stock
	double s_no_div_max = s_no_div0 * std::pow(u, N);
	for (int i = 0; i <= N; i++)
	{
		V[i] = std::max(-op.K + s_no_div_max, 0.0);
		s_no_div_max = s_no_div_max / u / u;
	}

	for (int k = N - 1; k >= 0; k--)
	{
		for (int i = 0; i <= k; i++)
		{
			V[i] = disc * (pu * V[i] + pd * V[i + 1]);
			if (k <= 2) option_result_vec[k][i] = V[i];
		}
	}

	// To calculate Greeks, we need to compute the first six corresponding stock prices at time 0, dt and 2 * dt
	// Integer coordinates corresponding to those dividend payment dates =  (dividend payment dates / dt)
	std::vector<double> div_date_vec = { 0 };
	for (int i = 0; i < num_div; i++)
	{
		int tmp = op.div_vec[i][2] / dt;
		div_date_vec.push_back(tmp);
	}

	// Let's say a stock pays discrete dividends at 0 < t1 < t2 < T, then the time between 0 and t1 is period 1
	auto get_present_period = [&](int i) {for (int j = num_div; j >= 0; j--) { if (i > div_date_vec[j]) return j; } return 0; };
	for (int k = 2; k >= 0; k--)
	{
		double s_no_div_max_ = s_no_div0 * std::pow(u, k);
		int present_period = get_present_period(k);
		for (int i = 0; i <= k; i++)
		{
			s_div_vec[i] = s_no_div_max_;
			for (int j = num_div - 1; j >= present_period; j--)
			{
				if (op.div_vec[j][0] == 0) { s_div_vec[i] /= (1 - op.div_vec[j][1]); }
				else { s_div_vec[i] += op.div_vec[j][1] * std::exp(-op.r * (op.div_vec[j][2] - k * dt)); }
			}
			stock_result_vec[k][i] = s_div_vec[i];
			s_no_div_max_ = s_no_div_max_ / u / u;
		}
	}

	// compute the Greeks
	double delta = (option_result_vec[1][0] - option_result_vec[1][1]) / (stock_result_vec[1][0] - stock_result_vec[1][1]);
	double gamma = ((option_result_vec[2][0] - option_result_vec[2][1]) / (stock_result_vec[2][0] - stock_result_vec[2][1]) -
		(option_result_vec[2][1] - option_result_vec[2][2]) / (stock_result_vec[2][1] - stock_result_vec[2][2])) / (stock_result_vec[2][0] - stock_result_vec[2][2]);
	double theta = (option_result_vec[2][1] - option_result_vec[0][0]) / (2 * dt);
	return { V[0], delta, gamma, theta };
}

// binomial tree pricer for european put option with discrete dividend payment with given accuracy
std::vector<double> european_call_binomial_tree(Option& op, int N0, double tol = 1e-4)
{
	double old_value = european_call_fixed_depth_binomial_tree(op, N0)[0];
	int N = N0 * 2;
	double new_value = european_call_fixed_depth_binomial_tree(op, N)[0];
	std::vector<double> tmp;
	while (std::abs(new_value - old_value) > tol)
	{
		old_value = new_value;
		N *= 2;
		tmp = european_call_fixed_depth_binomial_tree(op, N);
		new_value = tmp[0];
	}
	tmp.push_back(N);
	return tmp;
}
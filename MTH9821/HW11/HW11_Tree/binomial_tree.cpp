#include "binomial_tree.h"
#include<cmath>
#include<algorithm>

std::vector<double> BinomialTree(Option& op, int N0, double tol)
{
	if (op.type == European)
	{
		if (op.cp == call)
		{
			return european_call_binomial_tree(op, N0, tol);
		}
		else
		{
			return european_put_binomial_tree(op, N0, tol);
		}
	}
	else
	{
		if (op.cp == call)
		{
			return american_call_binomial_tree(op, N0, tol);
		}
		else
		{
			return american_put_binomial_tree(op, N0, tol);
		}
	}
}
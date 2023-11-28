#include<cstdlib>
#include<iostream>
#include<fstream>
#include<iomanip>
#include "binomial_tree.h"
int main()
{
	Option op;
	op.S0 = 50;
	op.K = 55.55;
	op.sigma = 0.25;
	op.r = 0.04;
	std::ofstream out;
    out.open("data.csv");

	// Question 1
	op.T = 3.0 / 12.0;
	op.div_vec = { {0, 0.015, 2.0 / 12.0} };

	auto tmp = european_put_binomial_tree(op, 6, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	tmp = american_put_binomial_tree(op, 6, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	op.div_vec = { {0, 0.015, 2.0 / 12.0}, {0, 0.015, 4.0 / 12.0}, {0, 0.015, 6.0 / 12.0} };
	op.T = 7.0 / 12.0;
	tmp = european_put_binomial_tree(op, 7, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	tmp = american_put_binomial_tree(op, 7, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;
	out << std::endl << std::endl;
	std::cout << "Question 1 Finished" << std::endl;

	// Question 2
	op.T = 3.0 / 12.0;
	op.div_vec = { {1, 0.95, 2.0 / 12.0} };
	tmp = european_put_binomial_tree(op, 6, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	tmp = american_put_binomial_tree(op, 6, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	op.div_vec = { {1, 0.95, 2.0 / 12.0}, {1, 0.95, 4.0 / 12.0}, {1, 0.95, 6.0 / 12.0} };
	op.T = 7.0 / 12.0;
	tmp = european_put_binomial_tree(op, 7, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	tmp = american_put_binomial_tree(op, 7, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;
	out << std::endl << std::endl;
	std::cout << "Question 2 Finished" << std::endl;

	// Question 3
	op.div_vec = { {1, 0.75, 2.0 / 12.0}, { 0, 0.02, 4.0 / 12.0 }, { 1, 0.25, 6.0 / 12.0 } };
	tmp = european_put_binomial_tree(op, 7, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;

	tmp = american_put_binomial_tree(op, 7, 1e-4);
	out << std::setprecision(10) << tmp[0] << "," << tmp[1] << "," << tmp[2] << "," << tmp[3] << "," << tmp[4] << std::endl;
	std::cout << "Question 3 Finished" << std::endl;
	out.close();
	system("pause");
	return 0;
}
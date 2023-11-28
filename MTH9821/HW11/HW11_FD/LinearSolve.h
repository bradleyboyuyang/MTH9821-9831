#pragma once
#include "Option.h"

class LinearSolver_Tridiag {
public:
	vector<vector<double>> L, U;

	LinearSolver_Tridiag(vector<vector<double>> A) {
		LU_tridiag(A);
	}

	void reset(vector<vector<double>> A) {
		LU_tridiag(A);
	}

	vector<double> solve(vector<double> b) {
		auto y = forward_solve(b);
		auto x = backward_solve(y);
		return x;
	}

private:
	void LU_tridiag(vector<vector<double>>& A) {
		int n = A.size();
		L.resize(n);
		for (auto & v : L) {
			v.resize(n, 0);
		}
		U.resize(n);
		for (auto & v : U) {
			v.resize(n, 0);
		}

		for (int i = 0; i < n - 1; ++i) {
			U[i][i] = A[i][i];
			U[i][i + 1] = A[i][i + 1];
			L[i][i] = 1;
			L[i + 1][i] = A[i + 1][i] / U[i][i];
			A[i + 1][i + 1] -= L[i + 1][i] * U[i][i + 1];
		}
		L[n-1][n-1] = 1;
		U[n - 1][n - 1] = A[n - 1][n - 1];
		return;
	}

	vector<double> forward_solve(vector<double> b) {
		vector<double> y(b.size());
		y[0] = b[0];
		for (int i = 1; i < b.size(); ++i) {
			y[i] = b[i] - L[i][i - 1] * y[i - 1];
		}
		return y;
	}

	vector<double> backward_solve(vector<double> y) {
		vector<double> x(y.size());
		int n = y.size();
		x[n - 1] = y[n - 1] / U[n - 1][n - 1];
		for (int i = n - 2; i >= 0; --i) {
			x[i] = (y[i] - U[i][i + 1] * x[i + 1]) / U[i][i];
		}
		return x;
	}

};




class LinearTimer_Tridiag {
public:
	vector<vector<double>> B;

	LinearTimer_Tridiag(vector<vector<double>> B_):B(B_) {}

	vector<double> time(vector<double> b) {
		int n = b.size();
		vector<double> res(n, 0);
		res[0] = B[0][0] * b[0] + B[0][1] * b[1];
		for (int i = 1; i < n - 1; i++) {
			res[i] = B[i][i - 1] * b[i - 1] + B[i][i] * b[i] + B[i][i + 1] * b[i + 1];
		}
		res[n - 1] = B[n - 1][n - 2] * b[n - 2] + B[n - 1][n - 1] * b[n - 1];
		return res;
	}

};

double Vec_norm(vector<double> a) {
	double res = 0.0;
	for (auto item : a) {
		res += item * item;
	}
	res = sqrt(res);
	return res;
}

vector<double> Vec_minus(vector<double> a, vector<double> b) {
	vector<double> res;
	for (int i = 0; i < a.size(); ++i) {
		res.push_back(a[i] - b[i]);
	}
	return res;
}

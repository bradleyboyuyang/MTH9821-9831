#pragma once
#include "Option.h"

class Generator{
public:
	virtual double operator()() = 0;
	virtual int countNum() { return 0; }
};

class LCG:public Generator {
public:
	long long seed;
	long a;
	long c;
	long k;
	
	LCG(long seed_=1, long a_=39373, long c_=0, long k_= 2147483647):
		seed(seed_),a(a_),c(c_),k(k_){}
	
	double operator()() {
		seed = (a * seed + c) % k;          //x_(i+1) = (a * x_i + c) mod k
		return (double)seed / (double)k;    //u_(i+1) = x_(i+1) / k
	}

};

class InverseTransform :public Generator {
	double a0 = 2.50662823884, a1 = -18.61500062529, a2 = 41.39119773534, a3 = -25.44106049637;
	double b0 = -8.47351093090, b1 = 23.08336743743, b2 = -21.06224101826, b3 = 3.13082909833;
	double c0 = 0.3374754822726147, c1 = 0.9761690190917186, c2 = 0.1607979714918209,
		c3 = 0.0276438810333863, c4 = 0.0038405729373609, c5 = 0.0003951896511919,
		c6 = 0.0000321767881768, c7 = 0.0000002888167364, c8 = 0.0000003960315187;

public:
	int count;
	LCG U;
	InverseTransform(LCG U_ = LCG()) :U(U_) { count = 0; }

	void reset() {
		count = 0;
		U.seed = 0;
	}

	int countNum() { return count; }

	double operator()() {
		double u = U();
		double y = u - 0.5;
		double r, x;
		if (abs(y) < 0.42) {
			r = y * y;
			x = y * (((a3*r + a2)*r + a1)*r + a0) / ((((b3*r + b2)*r + b1)*r + b0)*r + 1);
		}
		else {
			r = u;
			if (y > 0)r = 1 - u;
			r = log(-log(r));
			x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
			if (y < 0)x = -x;
		}
		++count;
		return x;
	}

};

class AcceptReject :public Generator {
	double c;
public:
	int count;
	LCG U;
	AcceptReject(LCG U_ = LCG()) :U(U_) { c = sqrt(2 * E / PI);count = 0; }

	void reset() {
		count = 0;
		U.seed = 0;
	}

	int countNum() { return count; }

	double operator()(){
		double u1 = U();
		double u2 = U();
		double u3 = U();
		double x = -log(u1);
		count += 3;
		while (u2 > exp(-0.5*(x - 1)*(x - 1))) {
			u1 = U();u2 = U();u3 = U();count += 3;
			x = -log(u1);
		}
		if (u3 <= 0.5)x = -x;
		return x;
	}
};

class BoxMuller :public Generator {
	bool has_store;
	double store;
public:
	int count;
	LCG U;
	BoxMuller(LCG U_ = LCG()) :U(U_) { count = 0; has_store = false; }

	void reset() {
		count = 0;
		U.seed = 0;
	}

	void refesh() {
		U.seed = 1;
		return;
	}

	int countNum() { return count; }

	double operator()() {
		if (has_store) {
			has_store = false;
			return store;
		}
		double u1 = U();
		double u2 = U();
		count += 2;
		u1 = 2 * u1 - 1;
		u2 = 2 * u2 - 1;
		double x = u1 * u1 + u2 * u2;
		while (x > 1) {
			u1 = U();u2 = U();count += 2;
			u1 = 2 * u1 - 1;
			u2 = 2 * u2 - 1;
			x = u1 * u1 + u2 * u2;
		}
		double y = sqrt(-2 * log(x) / x);
		double z1 = u1 * y, z2 = u2 * y;
		has_store = true;
		store = z2;
		return z1;
	}
};

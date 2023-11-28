#include "FiniteDifferencesPricer.h"

class DividendFiniteDifferencesPricer {
public:
	Option op;
	double x_compute;
	double x_left;
	double x_right;

	double alpha_temp;
	double alpha;
	int M;
	int N;

	double delta_tau;
	double delta_x;
	double tau_final;

	double a;
	double b;

	vector<vector<double>> Domin;// domin[n][m], first is space, second is time
	vector<double> x_pos;
	vector<double> tau_pos;

	//vector<vector<bool>> early_ex;

	double V_exact;
	//double omega, tol;

	FDType fd_type;

	vector<double> dividend;
	double tau_div;
	double x_compute_bar;
	double N_left;

	vector<vector<double>> Domin1;
	vector<double> x1_pos;
	vector<double> tau1_pos;
	double delta_tau1;
	double x_left_new;
	double x_right_new;
	double M1;
	double alpha1;

	DividendFiniteDifferencesPricer (Option op_, double alpha_temp_, int M_, FDType fdtype_,
	                       //double value_=0.0,// value is to set exact value. Used for American option
	                       //double omega_=1.2, double tol_=1e-6// omega and tol for SOR
	                       vector<double> dividend_)//{ratio,time}
		:op(op_), M(M_), alpha_temp(alpha_temp_),fd_type(fdtype_),dividend(dividend_) {

		//if (op.type == European||(op.type==DaO && op.cp==call))
			V_exact = op.BSvalue;
		//else if (op.type == American)
		//	V_exact = value_;
		/*
				if (fdtype_ == CN_SOR) {
					omega = omega_;
					tol = tol_;
				}
		*/
		// Domin1 parameter
		double x_left_temp;
		double x_right_temp;
		double N_right;
		tau_div=(op.T-dividend[1])*op.sigma*op.sigma/2.0;
		tau_final = op.T*op.sigma*op.sigma / 2.0;
		x_compute_bar=log(op.S0 / op.K)+log(1-dividend[0]);

		delta_tau=tau_div/(double)M;
		delta_x=sqrt(delta_tau/alpha_temp);

		        double temp = log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T;
		x_left_temp=temp - 3.0*op.sigma*sqrt(op.T);
		x_right_temp=temp + 3.0*op.sigma*sqrt(op.T);

		N_left=ceil((x_compute_bar-x_left_temp)/delta_x);
		N_right=ceil((x_right_temp-x_compute_bar)/delta_x);
		N=N_left+N_right;

		x_left=x_compute_bar-N_left*delta_x;
		x_right=x_compute_bar+N_right*delta_x;

		x_compute = log(op.S0 / op.K);
		tau_final = op.T*op.sigma*op.sigma / 2.0;

		alpha = alpha_temp;

		// change variable a, b
		a = (op.r - op.q) / op.sigma / op.sigma - 0.5;
		b = ((op.r - op.q) / op.sigma / op.sigma + 0.5)*((op.r - op.q) / op.sigma / op.sigma
		        + 0.5) + 2.0*op.q / op.sigma / op.sigma;

		// Domin
		for (int i = 0; i <= N; ++i) {
			vector<double> v(M + 1);
			Domin.emplace_back(v);
		}

		double x_temp = x_left;
		for (int i = 0; i <= N; ++i) {
			x_pos.emplace_back(x_temp);
			x_temp += delta_x;
		}

		double tau_temp = 0;
		for (int i = 0; i <= M; ++i) {
			tau_pos.emplace_back(tau_temp);
			tau_temp += delta_tau;
		}

		// Domin2
		x_left_new=x_left-log(1-dividend[0]);
		x_right_new=x_right-log(1-dividend[0]);
		M1=ceil((tau_final-tau_div)/alpha/delta_x/delta_x);
		delta_tau1=(tau_final-tau_div)/M1;

		for (int i = 0; i <= N; ++i) {
			vector<double> v(M1 + 1);
			Domin1.emplace_back(v);
		}

		x_temp = x_left_new;
		for (int i = 0; i <= N; ++i) {
			x1_pos.emplace_back(x_temp);
			x_temp += delta_x;
		}

		tau_temp -= delta_tau;
		for (int i = 0; i <= M1; ++i) {
			tau1_pos.emplace_back(tau_temp);
			tau_temp += delta_tau1;
		}
		alpha1=delta_tau1/delta_x/delta_x;

		// compute
		computeDomin();

		// measurements
		//computeVal();
		//computeRMS()
		computeGreeks();
		//computeRed();
		//computeEarly();

		return;
	}




	// main function
	void computeDomin() {
		// initialize
		for (int n = 0; n <= N; ++n) {
			Domin[n][0] = f(x_pos[n]);
		}
		for (int m = 1; m <= M; ++m) {
			Domin[0][m] = gleft(tau_pos[m]);
			Domin[N][m] = gright(tau_pos[m]);
		}


		// update rest domin
		if (fd_type == Forward) {
			for (int m = 0; m < M; ++m) {
				for (int n = 1; n < N; ++n) {
					if (op.type == European) {
						Domin[n][m + 1] = alpha * Domin[n - 1][m] + (1 - 2.0*alpha)*Domin[n][m]
						                  + alpha * Domin[n + 1][m];
					}
				}
			}
		}


		else if (fd_type == CN ) {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N - 1, 0.0));
			A[0][0] = 1 + alpha;
			A[0][1] = -alpha/2.0;
			for (int i = 1; i < N - 2; ++i) {
				A[i][i - 1] = -alpha/2.0;
				A[i][i] = 1 + alpha;
				A[i][i + 1] = -alpha/2.0;
			}
			A[N - 2][N - 3] = -alpha/2.0;
			A[N - 2][N - 2] = 1 + alpha;

			// initial B
			vector<vector<double>> B(N - 1, vector<double>(N - 1, 0));
			B[0][0] = 1 - alpha;
			B[0][1] = alpha / 2.0;
			for (int i = 1; i < N - 2; ++i) {
				B[i][i - 1] = alpha / 2.0;
				B[i][i] = 1 - alpha;
				B[i][i + 1] = alpha / 2.0;
			}
			B[N - 2][N - 3] = alpha / 2.0;
			B[N - 2][N - 2] = 1 - alpha;

			// compute
			LinearSolver_Tridiag ls_A(A);
			LinearTimer_Tridiag lt_B(B);
			vector<double> u_new = Get_time_slice(0,0);
			for (int m = 0; m < M; m++) {
				auto u_old = lt_B.time(u_new);
				u_old[0] += alpha / 2.0*Domin[0][m+1] + alpha / 2.0*Domin[0][m];
				u_old[N - 2] += alpha / 2.0*Domin[N][m+1] + alpha / 2.0*Domin[N][m];
				u_new = ls_A.solve(u_old);
				Set_time_slice(u_new, m + 1,0);
			}
		}

//////////////////////////////////
//////////////////////////////////

		double temp=pow((1-dividend[0]),-a);
		for (int n = 0; n <= N; ++n) {
			Domin1[n][0] = Domin[n][M]*temp;
		}
		for(int m=1; m<=M1; ++m) {
			Domin1[0][m]=gleft(tau1_pos[m])*temp;
			Domin1[N][m]=gright(tau1_pos[m])*temp;
		}

		// update rest domin
		if (fd_type == Forward) {
			for (int m = 0; m < M1; ++m) {
				for (int n = 1; n < N; ++n) {
					if (op.type == European) {
						Domin1[n][m + 1] = alpha1 * Domin1[n - 1][m] + (1 - 2.0*alpha1)*Domin1[n][m]
						                   + alpha1 * Domin1[n + 1][m];
					}
				}
			}
		}


		else if (fd_type == CN ) {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N - 1, 0.0));
			A[0][0] = 1 + alpha1;
			A[0][1] = -alpha1/2.0;
			for (int i = 1; i < N - 2; ++i) {
				A[i][i - 1] = -alpha1/2.0;
				A[i][i] = 1 + alpha1;
				A[i][i + 1] = -alpha1/2.0;
			}
			A[N - 2][N - 3] = -alpha1/2.0;
			A[N - 2][N - 2] = 1 + alpha1;

			// initial B
			vector<vector<double>> B(N - 1, vector<double>(N - 1, 0));
			B[0][0] = 1 - alpha1;
			B[0][1] = alpha1 / 2.0;
			for (int i = 1; i < N - 2; ++i) {
				B[i][i - 1] = alpha1 / 2.0;
				B[i][i] = 1 - alpha1;
				B[i][i + 1] = alpha1 / 2.0;
			}
			B[N - 2][N - 3] = alpha1 / 2.0;
			B[N - 2][N - 2] = 1 - alpha1;

			// compute
			LinearSolver_Tridiag ls_A(A);
			LinearTimer_Tridiag lt_B(B);
			vector<double> u_new = Get_time_slice(0,1);
			for (int m = 0; m < M1; m++) {
				auto u_old = lt_B.time(u_new);
				u_old[0] += alpha1 / 2.0*Domin1[0][m+1] + alpha1 / 2.0*Domin1[0][m];
				u_old[N - 2] += alpha1 / 2.0*Domin1[N][m+1] + alpha1 / 2.0*Domin1[N][m];
				u_new = ls_A.solve(u_old);
				Set_time_slice(u_new, m + 1,1);
			}
		}
		return;
	}


	/*
	 * measurements
	 */
/*
	// pointwise errors
	double v_approx;
	double v_approx2;

	double error_pointwise;
	double error_pointwise2;

	void computeVal() {
		int i = int((x_compute - x_left) / delta_x);
		double x_i = x_pos[i];
		double x_ip = x_pos[i + 1];
		double V_i = exp(-a * x_i - b * tau_final)*Domin[i][M];
		double V_ip = exp(-a * x_ip - b * tau_final)*Domin[i+1][M];

		double S_i = op.K*exp(x_i);
		double S_ip = op.K*exp(x_ip);

		// error1
		v_approx = ((S_ip - op.S0)*V_i + (op.S0 - S_i)*V_ip) / (S_ip - S_i);
		error_pointwise = abs(v_approx - V_exact);

		// error2
		double u_approx = ((x_ip - x_compute)*Domin[i][M] + (x_compute - x_i)*Domin[i + 1][M])
		                  / (x_ip - x_i);
		v_approx2 = exp(-a * x_compute - b * tau_final)*u_approx;
		error_pointwise2 = abs(v_approx2 - V_exact);

		return;
	}
*/
/*
	// RMS: root mean square error
	double RMS;

	void computeRMS() {
		double error_total = 0.0;
		int N_RMS = 0;
		for (int i = 0; i <= N; ++i) {
			double Sk = op.K*exp(x_pos[i]);
			double Vk = exp(-a * x_pos[i] - b * tau_final)*Domin[i][M];
			double Vexact = Option(op.type, op.cp, Sk, op.K, op.T, op.sigma, op.q, op.r)
			                .BSvalue;
			if (Vexact > 0.00001*op.S0) {
				N_RMS += 1;
				error_total += (Vk - Vexact)*(Vk - Vexact) / Vexact / Vexact;
			}
		}
		RMS = sqrt(error_total / (double)N_RMS);
		return;
	}

*/
	// Greeks error
	double v_approx;
	double delta_approx;
	double gamma_approx;
	double theta_approx;
	double u_value;

	void computeGreeks() {
		u_value=Domin1[N_left][M1];
		v_approx=exp(-a*x_compute-b*tau_final)*Domin1[N_left][M1];

		double x_i = x_compute;
		double x_ip = x_compute+delta_x;
		double x_im = x_compute-delta_x;
		//double x_ipp = x_compute+delta_x*2;

		double V_i = exp(-a * x_i - b * tau_final)*Domin1[N_left][M1];
		double V_ip = exp(-a * x_ip - b * tau_final)*Domin1[N_left+1][M1];
		double V_im = exp(-a * x_im - b * tau_final)*Domin1[N_left-1][M1];
		//double V_ipp = exp(-a * x_ipp - b * tau_final)*Domin[N_left+2][M1];

		double S_i = op.K*exp(x_i);
		double S_ip = op.K*exp(x_ip);
		double S_im = op.K*exp(x_im);
		//double S_ipp = op.K*exp(x_ipp);

		delta_approx = (V_ip - V_im) / (S_ip - S_im);
		gamma_approx = ((S_i-S_im)*V_ip-(S_ip-S_im)*V_i+(S_ip-S_i)*V_im)/(S_i-S_im)/(S_ip-S_i)/((S_ip-S_im)/2.0);
		

		double V_i_dt= exp(-a * x_i - b * (tau_final-delta_tau1))*Domin1[N_left][M1-1];
		//double V_ip_dt= exp(-a * x_ip - b * (tau_final - delta_tau))*Domin[i + 1][M-1];

		//double v_approx_dt = ((S_ip - op.S0)*V_i_dt + (op.S0 - S_i)*V_ip_dt) / (S_ip - S_i);
		theta_approx = (-V_i_dt + v_approx) / (2.0*delta_tau1 / op.sigma / op.sigma);

		return;
	}

/*
	// variance reduction method
	double v_approx_red;
	double error_pointwise_red;

	void computeRed() {
		if (op.type == American) {
			Option optemp(European, op.cp, op.S0, op.K, op.T, op.sigma, op.q, op.r);
			double v_approx_euro = FiniteDifferencePircer(optemp, alpha_temp, M, fd_type)
			                       .v_approx;
			double v_bs = op.BSvalue;
			v_approx_red = v_approx + (v_bs - v_approx_euro);
			error_pointwise_red = abs(v_approx_red - V_exact);
		}
		return;
	}


	// early exercise domin
	vector<double> early_t; // t from T to 0,
	vector<double> early_s; // at time t, the stock price which is optimal to early exercise

	void computeEarly() {
		if (op.type == American) {
			for (int m = 0; m <= M; ++m) {
				double t = op.T - 2 * m*delta_tau / op.sigma / op.sigma;
				early_t.emplace_back(t);
				// find stock price which is the first to early exercise
				int n;
				for (n = N; n > 1; --n) {
					if (early_ex[n - 1][m])
						break;
				}
				// n == Nopt+1
				double S_Nopt = op.K*exp(x_pos[n - 1]);
				double S_Nopt_p = op.K*exp(x_pos[n]);
				early_s.emplace_back(0.5*(S_Nopt + S_Nopt_p));
			}
		}
		return;
	}
*/

private:
	/*
	 * class base functions: do not modify
	 */

	// get vector of (U_1,U_2,...,U_N-1), time period m
	vector<double> Get_time_slice(int m, int j) {
		vector<double> res(N - 1);
		if(j==0) {
			for (int i = 1; i < N; ++i)
				res[i-1] = Domin[i][m];
		}
		if(j==1) {
			for (int i = 1; i < N; ++i)
				res[i-1] = Domin1[i][m];
		}

		return res;
	}

	void Set_time_slice(vector<double> u, int m,int j) {
		if(j==0) {
			for (int i = 1; i < N; ++i)
				Domin[i][m] = u[i-1];
		}
		if(j==1) {
			for (int i = 1; i < N; ++i)
				Domin1[i][m] = u[i-1];
		}
		return;
	}


	/*
	 * Bonudary functions: you can modify
	 */
	// Bonudary functions
	double f(double x) {
		if(op.cp==put)
			return op.K*exp(a*x)*max(1.0 - exp(x), 0.0);
		else if(op.cp==call)
			return op.K*exp(a*x)*max(exp(x)-1.0, 0.0);
	}

	double gleft(double tau) {
		if(op.cp==put) {
			if (op.type == European)
				return op.K*exp(a*x_left + b * tau)*(exp(-2.0*op.r*tau / op.sigma / op.sigma) -
				                                     exp(x_left - 2.0*op.q*tau / op.sigma / op.sigma));
			else if (op.type == American)
				return op.K*exp(a*x_left + b * tau)*(1.0 - exp(x_left));
		}
		if(op.cp==call) {
			return 0.0;
		}

	}

	double gright(double tau) {
		if(op.cp==put) {
			return 0.0;
		}
		if(op.cp==call ) {
			return op.K*exp(a*x_right + b * tau)*(-exp(-2.0*op.r*tau / op.sigma / op.sigma) +
			                                      exp(x_right - 2.0*op.q*tau / op.sigma / op.sigma));
		}

	}



	/*
	 * Domin computing functions: you can modify as you like
	 */
	/*
		void computeDomin_Forward() {
			for (int m = 0; m < M; ++m) {
				for (int n = 1; n < N; ++n) {
					if (op.type == European) {
						Domin[n][m + 1] = alpha * Domin[n - 1][m] + (1 - 2.0*alpha)*Domin[n][m]
						                  + alpha * Domin[n + 1][m];
					} else if (op.type == American) {
						double early_ex_premium = op.K*exp(a*x_pos[n] + b * tau_pos[m + 1])
						                          *max(1 - exp(x_pos[n]), 0.0);

						Domin[n][m + 1] = alpha * Domin[n - 1][m] + (1 - 2.0*alpha)*Domin[n][m]
						                  + alpha * Domin[n + 1][m];

						if (early_ex_premium > Domin[n][m + 1]) {
							early_ex[n][m + 1] = true;
							Domin[n][m + 1] = early_ex_premium;
						}
					}
				}
			}
			return;
		}
	*/
	/*
		void computeDomin_Backward() {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N-1,0.0));
			A[0][0] = 1 + 2 * alpha;
			A[0][1] = -alpha;
			for (int i = 1; i < N - 2; ++i) {
				A[i][i - 1] = -alpha;
				A[i][i] = 1 + 2 * alpha;
				A[i][i + 1] = -alpha;
			}
			A[N - 2][N - 3] = -alpha;
			A[N - 2][N - 2] = 1 + 2 * alpha;

			// compute
			LinearSolver_Tridiag ls(A);
			vector<double> u_new = Get_time_slice(0);
			for (int m = 0; m < M; m++) {
				auto u_old = u_new;
				u_old[0] += alpha * Domin[0][m + 1];
				u_old[N - 2] += alpha * Domin[N][m + 1];
				u_new = ls.solve(u_old);
				Set_time_slice(u_new, m + 1);
			}
			return;
		}
	*/
	/*
		void computeDomin_CN() {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N - 1, 0.0));
			A[0][0] = 1 + alpha;
			A[0][1] = -alpha/2.0;
			for (int i = 1; i < N - 2; ++i) {
				A[i][i - 1] = -alpha/2.0;
				A[i][i] = 1 + alpha;
				A[i][i + 1] = -alpha/2.0;
			}
			A[N - 2][N - 3] = -alpha/2.0;
			A[N - 2][N - 2] = 1 + alpha;

			// initial B
			vector<vector<double>> B(N - 1, vector<double>(N - 1, 0));
			B[0][0] = 1 - alpha;
			B[0][1] = alpha / 2.0;
			for (int i = 1; i < N - 2; ++i) {
				B[i][i - 1] = alpha / 2.0;
				B[i][i] = 1 - alpha;
				B[i][i + 1] = alpha / 2.0;
			}
			B[N - 2][N - 3] = alpha / 2.0;
			B[N - 2][N - 2] = 1 - alpha;

			// compute
			LinearSolver_Tridiag ls_A(A);
			LinearTimer_Tridiag lt_B(B);
			vector<double> u_new = Get_time_slice(0);
			for (int m = 0; m < M; m++) {
				auto u_old = lt_B.time(u_new);
				u_old[0] += alpha / 2.0*Domin[0][m+1] + alpha / 2.0*Domin[0][m];
				u_old[N - 2] += alpha / 2.0*Domin[N][m+1] + alpha / 2.0*Domin[N][m];
				u_new = ls_A.solve(u_old);
				Set_time_slice(u_new, m + 1);
			}
			return;
		}
	*/
	/*
		void computeDomin_CN_SOR() {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N - 1, 0.0));
			A[0][0] = 1 + alpha;
			A[0][1] = -alpha/2.0;
			for (int i = 1; i < N - 2; ++i) {
				A[i][i - 1] = -alpha/2.0;
				A[i][i] = 1 + alpha;
				A[i][i + 1] = -alpha/2.0;
			}
			A[N - 2][N - 3] = -alpha/2.0;
			A[N - 2][N - 2] = 1 + alpha;

			// initial B
			vector<vector<double>> B(N - 1, vector<double>(N - 1, 0));
			B[0][0] = 1 - alpha;
			B[0][1] = alpha / 2.0;
			for (int i = 1; i < N - 2; ++i) {
				B[i][i - 1] = alpha / 2.0;
				B[i][i] = 1 - alpha;
				B[i][i + 1] = alpha / 2.0;
			}
			B[N - 2][N - 3] = alpha / 2.0;
			B[N - 2][N - 2] = 1 - alpha;

			// compute
			LinearTimer_Tridiag lt_B(B);

			for (int m = 0; m < M; m++) {
				vector<double> u = Get_time_slice(m);
				u = lt_B.time(u);
				u[0] += alpha / 2.0*Domin[0][m + 1] + alpha / 2.0*Domin[0][m];
				u[N - 2] += alpha / 2.0*Domin[N][m + 1] + alpha / 2.0*Domin[N][m];
				// u now is b


				vector<double> early_g;
				// calculate early premium for American Options
				if (op.type == American) {
					for (int i = 1; i < N; i++) {
						early_g.push_back(op.K*exp(a*x_pos[i] + b * tau_pos[m + 1])
						                  *max(1 - exp(x_pos[i]), 0.0));
					}
				}


				vector<double> x_old(N - 1, -1.0);
				vector<double> x_new(N-1,0.0);
				if(op.type==American)x_new=early_g;//initial guess
				while (Vec_norm(Vec_minus(x_old, x_new)) >= tol) {
					x_old = x_new;
					for (int i = 0; i<N - 1; ++i) {
						double x1 = 0, x2 = 0;
						if (i > 0)x1 = x_new[i - 1];
						if (i < N - 2)x2 = x_old[i + 1];
						x_new[i] = (1 - omega)*x_old[i] + omega * alpha / 2.0 / (1 + alpha)*(x2+x1)
						           + omega / (1 + alpha)*u[i];
						// if American
						if(op.type==American)
							x_new[i] = max(x_new[i], early_g[i]);
					}
				}

				Set_time_slice(x_new, m + 1);
			}
			return;
		}

	*/

};

#pragma once
#include"Option.h"
#include"LinearSolve.h"
#include <math.h>

enum FDType { Forward, Backward, CN, CN_SOR };

// Only for Put!
class FiniteDifferencePircer {
public:
	DAO_Call_Euro_Option op;
    
    
	double x_compute;
	double x_left;
	double x_right;
    double x_right_tilta;

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

	vector<vector<bool>> early_ex;
	
	double V_exact;
	double omega, tol;

	FDType fd_type;
    
    double delta_x_temp;
    int N_left;
    int N_right;

	FiniteDifferencePircer(DAO_Call_Euro_Option op_, double alpha_temp_, int M_, FDType fdtype_,
		double value_=0.0,// value is to set exact value. Used for American option
		double omega_=1.2, double tol_=1e-6)// omega and tol for SOR
		:op(op_), M(M_), alpha_temp(alpha_temp_),fd_type(fdtype_) {

            if (op.type == European)
                V_exact = op.BSvalue;
            else if (op.type == American)
                V_exact = value_;

            if (fdtype_ == CN_SOR) {
                omega = omega_;tol = tol_;
            }
            
            // Domin parameter
            x_compute = log(op.S0 / op.K);
            tau_final = op.T*op.sigma*op.sigma / 2.0;
            double temp = log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T;
            
            x_left = log(op.B / op.K);
            x_right_tilta = temp + 3.0*op.sigma*sqrt(op.T);
            
            //Computational Domain
            delta_tau = tau_final / (double)M;
            delta_x_temp = sqrt(delta_tau / alpha_temp);
            N_left = int((x_compute - x_left)/delta_x_temp);
            delta_x = (x_compute - x_left) / (double)N_left;
            alpha = delta_tau / delta_x / delta_x;
            N_right = ceil((x_right_tilta - x_compute)/delta_x);
            N = N_left + N_right;
            x_right = x_compute + N_right*delta_x;

            // change variable a, b
            a = (op.r - op.q) / op.sigma / op.sigma - 0.5;
            b = ((op.r - op.q) / op.sigma / op.sigma + 0.5)*((op.r - op.q) / op.sigma / op.sigma + 0.5) + 2.0*op.q / op.sigma / op.sigma;

            // Domin
            for (int i = 0;i <= N;++i) {
                vector<double> v(M + 1);
                Domin.emplace_back(v);
                vector<bool> w(M + 1, false);
                early_ex.emplace_back(w);
            }

            double x_temp = x_left;
            for (int i = 0;i <= N;++i) {
                x_pos.emplace_back(x_temp);
                x_temp += delta_x;
            }

            double tau_temp = 0;
            for (int i = 0;i <= M;++i) {
                tau_pos.emplace_back(tau_temp);
                tau_temp += delta_tau;
            }
            
            // compute
            computeDomin();
            computeVal();
            computeGreeks();
            return;
	}

	// Bonudary functions
	double f(double x) {
		return op.K*exp(a*x)*max(exp(x)-1, 0.0);
	}
    
    double gleft(double tau) {
        return 0.0;
    }
    
	double gright(double tau) {
		if (op.type == European)
			return op.K*exp(a*x_right + b * tau)*(-exp(-2.0*op.r*tau / op.sigma / op.sigma) +
				exp(x_right - 2.0*op.q*tau / op.sigma / op.sigma));
        
        return -1;
	}




	// main function
	void computeDomin() {
		// initialize
		for (int n = 0;n <= N;++n) {
			Domin[n][0] = f(x_pos[n]);
		}
		for (int m = 1;m <= M;++m) {
			Domin[0][m] = gleft(tau_pos[m]);
			Domin[N][m] = gright(tau_pos[m]);
		}

		// update rest domin
		if (fd_type == Forward) {
			for (int m = 0;m < M;++m) {
				for (int n = 1;n < N;++n) {
					if (op.type == European) {
						Domin[n][m + 1] = alpha * Domin[n - 1][m] + (1 - 2.0*alpha)*Domin[n][m]
							+ alpha * Domin[n + 1][m];
					}
					else if (op.type == American) {
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
		}

		else if (fd_type == Backward) {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N-1,0.0));
			A[0][0] = 1 + 2 * alpha;A[0][1] = -alpha;
			for (int i = 1;i < N - 2;++i) {
				A[i][i - 1] = -alpha;A[i][i] = 1 + 2 * alpha;A[i][i + 1] = -alpha;
			}
			A[N - 2][N - 3] = -alpha;A[N - 2][N - 2] = 1 + 2 * alpha;

			// compute
			LinearSolver_Tridiag ls(A);
			vector<double> u_new = Get_time_slice(0);
			for (int m = 0;m < M;m++) {
				auto u_old = u_new;
				u_old[0] += alpha * Domin[0][m + 1];
				u_old[N - 2] += alpha * Domin[N][m + 1];
				u_new = ls.solve(u_old);
				Set_time_slice(u_new, m + 1);
			}
		}


		else if (fd_type == CN || fd_type==CN_SOR) {
			// initial A
			vector<vector<double>> A(N - 1, vector<double>(N - 1, 0.0));
			A[0][0] = 1 + alpha;A[0][1] = -alpha/2.0;
			for (int i = 1;i < N - 2;++i) {
				A[i][i - 1] = -alpha/2.0;A[i][i] = 1 + alpha;A[i][i + 1] = -alpha/2.0;
			}
			A[N - 2][N - 3] = -alpha/2.0;A[N - 2][N - 2] = 1 + alpha;

			// initial B
			vector<vector<double>> B(N - 1, vector<double>(N - 1, 0));
			B[0][0] = 1 - alpha;B[0][1] = alpha / 2.0;
			for (int i = 1;i < N - 2;++i) {
				B[i][i - 1] = alpha / 2.0;B[i][i] = 1 - alpha;B[i][i + 1] = alpha / 2.0;
			}
			B[N - 2][N - 3] = alpha / 2.0;B[N - 2][N - 2] = 1 - alpha;

			// compute
			// if CN
			if (fd_type == CN) {
				LinearSolver_Tridiag ls_A(A);
				LinearTimer_Tridiag lt_B(B);
				vector<double> u_new = Get_time_slice(0);
				for (int m = 0;m < M;m++) {
					auto u_old = lt_B.time(u_new);
					u_old[0] += alpha / 2.0*Domin[0][m+1] + alpha / 2.0*Domin[0][m];
					u_old[N - 2] += alpha / 2.0*Domin[N][m+1] + alpha / 2.0*Domin[N][m];
					u_new = ls_A.solve(u_old);
					Set_time_slice(u_new, m + 1);
				}
			}

			// if CN_SOR
			else if (fd_type == CN_SOR) {
				LinearTimer_Tridiag lt_B(B);
			
				for (int m = 0;m < M;m++) {
					vector<double> u = Get_time_slice(m);
					u = lt_B.time(u);
					u[0] += alpha / 2.0*Domin[0][m + 1] + alpha / 2.0*Domin[0][m];
					u[N - 2] += alpha / 2.0*Domin[N][m + 1] + alpha / 2.0*Domin[N][m];
					// u now is b

					//European
					if (op.type == European) {
						vector<double> x_old(N - 1, 1.0);
						vector<double> x_new(N - 1, 0.0);//initial guess
						while (Vec_norm(Vec_minus(x_old, x_new))>=tol) {
							x_old = x_new;
							x_new[0]= (1 - omega)*x_old[0] + omega * alpha / 2.0 / (1 + alpha)*(x_old[1])
								+ omega / (1 + alpha)*u[0];
							for (int i = 1;i < N - 2;++i) {
								x_new[i] = (1 - omega)*x_old[i] + omega * alpha / 2.0 / (1 + alpha)*(x_old[i + 1]
									+ x_new[i - 1]) + omega / (1 + alpha)*u[i];
							}
							x_new[N-2] = (1 - omega)*x_old[N-2] + omega * alpha / 2.0 / (1 + alpha)*(x_new[N-3]) 
								+ omega / (1 + alpha)*u[N-2];
						}

						Set_time_slice(x_new, m + 1);
					}

					//American
					if (op.type == American) {
						vector<double> early_g;
						for (int i = 1;i < N;i++) {
							early_g.push_back(op.K*exp(a*x_pos[i] + b * tau_pos[m + 1])
								*max(1 - exp(x_pos[i]), 0.0));
						}

						vector<double> x_old(N - 1, -1.0);
						vector<double> x_new=early_g;//initial guess
						while (Vec_norm(Vec_minus(x_old, x_new)) >= tol) {
							x_old = x_new;
							for (int i = 0;i<N - 1;++i) {
								double x1 = 0, x2 = 0;
								if (i > 0)x1 = x_new[i - 1];
								if (i < N - 2)x2 = x_old[i + 1];
								x_new[i] = (1 - omega)*x_old[i] + omega * alpha / 2.0 / (1 + alpha)*(x2+x1)
									+ omega / (1 + alpha)*u[i];
								x_new[i] = max(x_new[i], early_g[i]);
							}
						}

						Set_time_slice(x_new, m + 1);
					}
					
				}
			}
			
		}

		

		return;
	}


	// get vector of (U_1,U_2,...,U_N-1), time period m
	vector<double> Get_time_slice(int m) {
		vector<double> res(N - 1);
		for (int i = 1;i < N;++i) {
			res[i-1] = Domin[i][m];
		}
		return res;
	}

	void Set_time_slice(vector<double> u, int m) {
		for (int i = 1;i < N;++i)
			Domin[i][m] = u[i-1];
		return;
	}

	double v_approx;
	//double v_approx2;

	double error_pointwise;
	//double error_pointwise2;
    
    double u_value;
    
	void computeVal() {
		int i = N_left;
        u_value = Domin[i][M];
		double V_i = exp(-a * x_compute - b * tau_final) * u_value;
		v_approx = V_i;
		error_pointwise = abs(v_approx - V_exact);

        /*
		double u_approx = ((x_ip - x_compute)*Domin[i][M] + (x_compute - x_i)*Domin[i + 1][M])
			/ (x_ip - x_i);
		v_approx2 = exp(-a * x_compute - b * tau_final)*u_approx;
		error_pointwise2 = abs(v_approx2 - V_exact);*/
		return;
	}


/*
	double RMS;

	void computeRMS() {
		double error_total = 0.0;
		int N_RMS = 0;
		for (int i = 0;i <= N;++i) {
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

	double delta_approx;
	double gamma_approx;
	double theta_approx;

	void computeGreeks() {
        V_exact = op.BSvalue;
        int i = N_left;
        double x_i_prev = x_pos[i - 1];
        double x_i = x_pos[i];
        double x_i_next = x_pos[i + 1];
        
        double V_i_prev = exp(-a * x_i_prev - b * tau_final) * Domin[i-1][M];
        double V_i = exp(-a * x_i - b * tau_final) * Domin[i][M];
        double V_i_next = exp(-a * x_i_next - b * tau_final) * Domin[i+1][M];
        
        double S_i_prev = op.K * exp(x_i_prev);
        double S_i = op.K * exp(x_i);
        double S_i_next = op.K * exp(x_i_next);
        
		delta_approx = (V_i_next - V_i_prev) / (S_i_next - S_i_prev);
		gamma_approx = ((S_i-S_i_prev)*V_i_next - (S_i_next-S_i_prev)*V_i + (S_i_next - S_i)*V_i_prev)/(S_i-S_i_prev)/(S_i_next-S_i)/((S_i_next-S_i_prev)/2);
        
        double dt = 2 * delta_tau/op.sigma/op.sigma;
		double V_approx_dt= exp(-a * x_i - b * (tau_final-delta_tau))*Domin[i][M-1];
        theta_approx = (v_approx - V_approx_dt) / dt;
		
		return;
	}
    
    
/*
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
*/

    /*
	vector<double> early_t;
	vector<double> early_s;

	void computeEarly() {
		if (op.type == American) {
			for (int m = 0;m <= M;++m) {
				double t = op.T - 2 * m*delta_tau / op.sigma / op.sigma;
				early_t.emplace_back(t);
				int n;
				for (n = N;n > 1;--n) {
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
};

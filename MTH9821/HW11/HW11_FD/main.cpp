#include "FiniteDifferencesPricer.h"
#include "DividendFiniteDifferencePricer.h"

int main() {
	Option op(European,call,53, 50, 1, 0.25, 0, 0.04);
	vector<double> dividend{0.01,5.0/12.0};
	//FiniteDifferencePircer pricer(op, 0.45, 4, CN_SOR, exact_value);

	ofstream out;
	out.open("1_data.csv");
	out<<"new1"<<endl<<endl;
	
	for(int M=4;M<=64*4;M*=4){
		DividendFiniteDifferencesPricer pricer(op,0.4,M,Forward,dividend);
		out<<M<<","
			<<setprecision(10)<<pricer.M1<<","<<pricer.alpha1<<","<<pricer.N<<","
			<<pricer.x_left<<","<<pricer.x_right<<","<<pricer.x_left_new<<","
			<<pricer.x_right_new<<","<<pricer.tau_div<<","<<pricer.delta_tau<<","
			<<pricer.delta_tau1<<","<<pricer.delta_x<<endl;
	}
	
	out<<endl<<endl<<endl;
	for(int M=4;M<=64*4;M*=4){
		DividendFiniteDifferencesPricer pricer(op,4,M,CN,dividend);
		out<<M<<","
			<<setprecision(10)<<pricer.M1<<","<<pricer.alpha1<<","<<pricer.N<<","
			<<pricer.x_left<<","<<pricer.x_right<<","<<pricer.x_left_new<<","
			<<pricer.x_right_new<<","<<pricer.tau_div<<","<<pricer.delta_tau<<","
			<<pricer.delta_tau1<<","<<pricer.delta_x<<endl;
	}
	
	
	
	out<<endl<<endl<<endl;
	for(int M=4;M<=64*4;M*=4){
		DividendFiniteDifferencesPricer pricer(op,0.4,M,Forward,dividend);
		out<<M<<","
			<<setprecision(10)<<pricer.u_value<<","<<pricer.v_approx<<","
			<<pricer.delta_approx<<","<<pricer.gamma_approx<<","
			<<pricer.theta_approx<<endl;
	}
	
	out<<endl<<endl<<endl;
	for(int M=4;M<=64*4;M*=4){
		DividendFiniteDifferencesPricer pricer(op,4,M,CN,dividend);
		out<<M<<","
			<<setprecision(10)<<pricer.u_value<<","<<pricer.v_approx<<","
			<<pricer.delta_approx<<","<<pricer.gamma_approx<<","
			<<pricer.theta_approx<<endl;
	}
	
	out<<endl<<endl<<endl;
	
	
	DividendFiniteDifferencesPricer pricer(op,0.4,4,Forward,dividend);
	for (int m = 0;m <= pricer.M;++m) {
		for (int n = 0;n <= pricer.N;++n) {
			out << setprecision(10) << pricer.Domin[n][m] << ",";
		}
		out << endl;
	}
	out << endl << endl;
	for (int m = 0;m <= pricer.M1;++m) {
		for (int n = 0;n <= pricer.N;++n) {
			out << setprecision(10) << pricer.Domin1[n][m] << ",";
		}
		out << endl;
	}
	
	out.close();
	system("pause");
	return 0;
}
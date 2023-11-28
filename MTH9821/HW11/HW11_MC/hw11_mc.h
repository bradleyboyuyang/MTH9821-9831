#ifndef hw11_mc_h
#define hw11_mc_h

#include "Option.h"
#include "Generator.h"
#include <vector>

void hw_11(){
    double S0,K,T,sigma,r,q;
    double t1,t2,t3;
    long long N,n,m;
    vector<double> vec_bm;
    double q_1_fx = 0.95;
    double q_2 = 0.02;
    double q_3_fx = 0.45;
    S0 = 50;
    K = 55.55;
    t1 = 2.0/12;
    t2 = 4.0/12;
    t3 = 6.0/12;
    T = 7.0/12;
    sigma = 0.25;
    r = 0.04;
    q = 0;
    N = 4*10000*pow(2,8);
    m = 4;
    n = N/m;
    
    //no_dividend paying put;
    Option p(European, put, S0,K,T,sigma,q,r);
    double p_bs = p.BSvalue;
    double delta_bs = p.BSdelta;
    
    BoxMuller bm;
    for(long i = 0; i < N; i++){
        vec_bm.push_back(bm());
    }
    
    //without variance reduction
    double pi,pi_tilta,si1,si2,si3,si4,si_tilta;
    double p_hat, p_tilta_hat, p_cv_hat;
    double wi,b_hat;
    vector<double> vec_p_1,vec_delta_1,vec_p_cv,vec_delta_cv;
    
    for(int k = 0; k <= 8; k++){
        N = m * 10000 * pow(2,k);
        n = N/m;
        pi = 0;
        double pi_sum = 0; double pi_tilta_sum = 0; double wi_sum = 0;
        vector<double> vec_pi;
        vector<double> vec_pi_tilta;
        vector<double> vec_siT;
        vector<double> vec_siT_tilta;
		
        // generate S_i with dividend and without dividend (for control variate)
        for(long i = 0; i < n; i++){
            //European put with three discrete dividends
            si1 = S0 * exp((r-sigma*sigma/2)*t1 + sigma * sqrt(t1)*vec_bm[m*i])-q_1_fx;
            si2 = si1 * exp((r-sigma*sigma/2)*(t2-t1) + sigma * sqrt(t2-t1)*vec_bm[m*i+1])*(1-q_2);
            si3 = si2 * exp((r-sigma*sigma/2)*(t3-t2) + sigma * sqrt(t3-t2)*vec_bm[m*i+2])-q_3_fx;
            si4 = si3 * exp((r-sigma*sigma/2)*(T-t3) + sigma * sqrt(T-t3)*vec_bm[m*i+3]);
            pi = exp(-r*T) * max(K-si4,0.0);
            pi_sum += pi;
            vec_pi.push_back(pi);
            vec_siT.push_back(si4);
            
            //use the same option on a non-dividend paying asset as a control variate
            double temp = sqrt(t1)*vec_bm[m*i] + sqrt(t2-t1)*vec_bm[m*i+1] + sqrt(t3-t2)*vec_bm[m*i+2] + sqrt(T-t3)*vec_bm[m*i+3];
            si_tilta = S0 * exp((r-sigma*sigma/2)*T + sigma * temp);
            pi_tilta = exp(-r*T) * max(K-si_tilta, 0.0);
            pi_tilta_sum += pi_tilta;
            vec_pi_tilta.push_back(pi_tilta);
            vec_siT_tilta.push_back(si_tilta);
        }
        // calculate average to calculate b_hat
        p_hat = pi_sum / (double) n;
        vec_p_1.push_back(p_hat);
        p_tilta_hat = pi_tilta_sum / (double) n;
        
        std::cout<<p_hat<<", "<< p_tilta_hat<<", "<<p_bs<<endl;
        
        double temp1 = 0; double temp2 = 0;
        
        for(long i = 0; i<n;i++){
            temp1 += (vec_pi_tilta[i]-p_tilta_hat)*(vec_pi[i]-p_hat);
            temp2 += (vec_pi_tilta[i]-p_tilta_hat)*(vec_pi_tilta[i]-p_tilta_hat);
        }
        b_hat = temp1/temp2;
		
        // calculate control variate price
        for(long i = 0; i<n;i++){
            wi = vec_pi[i] - b_hat*(vec_pi_tilta[i] - p_bs);
            wi_sum += wi;
        }
        p_cv_hat = wi_sum/(double)n;
        vec_p_cv.push_back(p_cv_hat);
        
        
        //delta
        double deltai = 0;
        double delta_hat = 0;
        double deltai_tilta = 0;
        double delta_tilta_hat = 0;
        double delta_cv_hat = 0;
        vector<double> vec_deltai;
        vector<double> vec_deltai_tilta;
        
		// calculate delta with and without dividend
        for(long i = 0; i < n; i++){
            //with dividend
            deltai = 0;deltai_tilta = 0;
            if (K>vec_siT[i]){
                deltai = -exp(-r*T) * vec_siT_tilta[i] / (double)S0 * (1-q_2);
            }
            delta_hat += deltai;
            vec_deltai.push_back(deltai);
            
            //Use again, as control variate,
            //the Delta of the option on the non-dividend-paying version of the asset.
            if (K > vec_siT_tilta[i]){
                deltai_tilta = -exp(-r*T) * vec_siT_tilta[i] / (double)S0;
            }
            delta_tilta_hat += deltai_tilta;
            vec_deltai_tilta.push_back(deltai_tilta);
        }
        
        delta_hat = delta_hat/(double)n;
        vec_delta_1.push_back(delta_hat);
        delta_tilta_hat = delta_tilta_hat/(double)n;
        
        
        // calculate b_hat for delta control variate
        temp1 = 0; temp2 = 0;wi_sum = 0;
        for(long i = 0; i<n;i++){
            temp1 += (vec_deltai_tilta[i]-delta_tilta_hat)*(vec_deltai[i]-delta_hat);
            temp2 += (vec_deltai_tilta[i]-delta_tilta_hat)*(vec_deltai_tilta[i]-delta_tilta_hat);
        }
        b_hat = temp1/temp2;
        
		// calculate cv delta
        for(long i = 0; i<n;i++){
            wi = vec_deltai[i] - b_hat*(vec_deltai_tilta[i] - delta_bs);
            wi_sum += wi;
        }
        delta_cv_hat = wi_sum/(double)n;
        vec_delta_cv.push_back(delta_cv_hat);
    }

    
    //output the result
    ofstream out;
    out.open("hw_11_2.csv");
    for(int k = 0; k<=8; k++){
        out<<setprecision(10)<<vec_p_1[k]<<", "<<vec_delta_1[k]<<endl;
    }
    out<<endl;
    
    
    for(int k = 0; k<=7; k++){
        out<<setprecision(10)<<vec_p_cv[k]<<", "<<vec_delta_cv[k]<<endl;
    }
    out<<endl;
    
    
    out.close();
    
}

#endif /* hw11_mc_h */

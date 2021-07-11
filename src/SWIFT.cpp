/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"SWIFT.h"
#include"BSM.h"
#include <iostream>
#include<fftw3.h>
#include<complex>

using namespace std::complex_literals;

#define TRUNCATION_PRECISION 0.0005

SWIFT::SWIFT(const HDistribution& distr, const ffloat stock_price, const options_chain& opts,const unsigned int m_d) : distr(distr), S(stock_price), m(m_d){
    ffloat tau=static_cast<ffloat>(opts.days_to_expiry)/trading_days;
    ffloat max=yearly_risk_free*tau+std::log(stock_price/opts.min_strike);
    ffloat min=yearly_risk_free*tau+std::log(stock_price/opts.max_strike);
    ffloat c=std::abs(distr.first_order_moment(tau))+10.*std::sqrt(std::fabs(distr.second_order_moment(tau))+std::sqrt(std::abs(distr.fourth_order_moment(tau))));
    this->exp2_m=std::exp2(m);
    this->lower=min-c;
    this->upper=max+c;
    this->k_1=ceil(exp2_m*(min - c));
    this->k_2=floor(exp2_m*(max + c));
    ffloat iota_density=ceil(std::log2(M_PI*std::abs(k_1-k_2)))-1;
    this->J=std::exp2(iota_density-1);
    std::cout<<"m: "<<m<<"\texp2_m: "<<exp2_m<<"\tc:"<<c<<"\tLower integral bound: "<<k_1<<"\tUpper integral bound: "<<k_2<<"\tJ: "<<J<<'\n';
    //int cas;
    //for(i=1;(i<=k_2-k_1)||(cas=(i<=-k_1));i++) payoffs[-k_1+cas*i-(1-cas)*(i+k_1))]=(out[i]+(cas*2-1)*out2[i-1])*0.5;
    
}
unsigned int SWIFT::get_m(const HDistribution& distr,ffloat tau){
    unsigned int i=0;
    while(distr.int_error(i++,tau)>TRUNCATION_PRECISION);
    return i;
}
std::complex<ffloat> SWIFT::H(ffloat y, int j){
    ffloat u_j=M_PI*(2*static_cast<ffloat>(j)+1)/(4.*J)*exp2_m;
    //std::cout<<"u_j"<<u_j<<'\n';
    return 1i*std::exp(1i*u_j*y)*(1/u_j-std::exp(y)/(-1i+u_j));
}
void SWIFT::get_FFT_coeffs(){
    int i,j;
    double *in,*out,*in2,*out2;
    fftw_plan plan1;
    fftw_plan plan2;
    
    double ea = 1.;
    double eb = exp(upper);
    
    double * d = (double*) fftw_malloc(sizeof(double)*J*2.);
    double * e = (double*) fftw_malloc(sizeof(double)*J*2.);
    
    double * real_input = (double*) fftw_malloc(sizeof(double)*J*2.);
    double * complex_input = (double*) fftw_malloc(sizeof(double)*J*2.);

    in = (double*) fftw_malloc(sizeof(double)*J*2.);
    out = (double*) fftw_malloc(sizeof(double)*J*2.);
    in2 = (double*) fftw_malloc(sizeof(double)*J*2.);
    out2 = (double*) fftw_malloc(sizeof(double)*J*2.);

    for(j=0;j<J;j++)
    {
        double Co=((2*j+1.)/(2.*J*2.))*M_PI;
        //std::cout<<"2u_j"<<Co<<'\n';
        double B=(1./(Co*exp2_m));
        double A=(Co*exp2_m)/(1.+pow(Co*exp2_m,2));
        double sb=std::sin(Co*exp2_m*upper);
        double sa=std::sin(Co*exp2_m*0.);
        double cb=std::cos(Co*exp2_m*upper);
        double ca=std::cos(Co*exp2_m*0.);
        double I11=eb*sb-ea*sa+B*eb*cb-B*ea*ca;
        double I12=-eb*cb+ea*ca+B*eb*sb-B*ea*sa;
        double I21=sb-sa;
        double I22=ca-cb;
        d[j]=A*I11-B*I21;
        e[j]=A*I12-B*I22;
        std::complex<ffloat> current=H(upper,j)-H(0,j);// pricing of call options
        real_input[j]=current.real();
        complex_input[j]=current.imag();
        //std::cout<<"eb"<<eb<<"real_input: "<<real_input[j]<<"\tcomplex_input: "<<complex_input[j]<<"\tERinput real: "<<d[j]<<"\tERinput compl: "<<e[j]<<"\tdiff real: "<<std::fabs(real_input[j]-d[j])<<"\tdiff compl: "<<std::fabs(complex_input[j]-e[j])<<'\n';
    }
    //std::cout<<"Start test\n";

    //Calculem amb FFT

    plan1 = fftw_plan_r2r_1d(2.*J, d, out, FFTW_REDFT10, FFTW_ESTIMATE);       //Here we set which kind of transformation we want to perform
    fftw_execute(plan1);                                                             //Execution of FFT

    plan2 = fftw_plan_r2r_1d(2.*J, e, out2, FFTW_RODFT10, FFTW_ESTIMATE);       //Here we set which kind of transformation we want to perform
    fftw_execute(plan2);                                                             //Execution of FFT


    std::vector<double> payoffs(k_2-k_1 + 1);

    //Para k=0
    payoffs[-k_1]=(1./J)*(out[0]/2.);
    //printf("res[%d]=%lf\n",-k_1,payoffs[-k_1]);

    //Para k>0
    
    for(i=1;i<=k_2;i++)
    {
        payoffs[-k_1+i]=(1./J)*((out[i]+out2[i-1])/2.); //out viene multiplicado por 2 !!!
        printf("res[%d]=%lf\n",-k_1+i,payoffs[-k_1+i]);
    }

    //Para k<0

    for(i=1;i<=-k_1;i++)
    {
        payoffs[-k_1-i]=(1./J)*((out[i]-out2[i-1])/2.); //out viene multiplicado por 2 !!!
        printf("res[%d]=%lf\n",-k_1-i,payoffs[-k_1-i]);
    }

    fftw_destroy_plan(plan1);                                            //Destroy plan
    fftw_destroy_plan(plan2);
    free(real_input);
    free(complex_input);
    free(d);
    free(e);
    free(in);                                                            //Free memory
    free(out);                                                           //Free memory
    free(in2);                                                           //Free memory
    free(out2);
}

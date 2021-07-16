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
    std::cout<<"min_strike: "<<opts.min_strike<<"\tmax_strike: "<<opts.max_strike<<'\n';
    ffloat max=yearly_risk_free*tau+std::log(stock_price/opts.min_strike);
    ffloat min=yearly_risk_free*tau+std::log(stock_price/opts.max_strike);
    ffloat c=std::abs(distr.first_order_moment(tau))+10.*std::sqrt(std::fabs(distr.second_order_moment(tau))+std::sqrt(std::abs(distr.fourth_order_moment(tau))));
    this->exp2_m=std::exp2(m);
    this->sqrt_exp2_m=std::sqrt(static_cast<ffloat>(exp2_m));
    this->lower=min-c;
    this->upper=max+c;
    this->k_1=ceil(exp2_m*lower);
    this->k_2=floor(exp2_m*upper);
    ffloat iota_density=ceil(std::log2(M_PI*std::abs(k_1-k_2)))-1;
    this->J=std::exp2(iota_density-1);
    density_coeffs=new std::vector<std::complex<ffloat>>(J);
    std::cout<<"m: "<<m<<"\texp2_m: "<<exp2_m<<"\tc:"<<c<<"\tLower integral bound: "<<k_1<<"\tUpper integral bound: "<<k_2<<"\tJ: "<<J<<'\n';
    this->get_FFT_coeffs();
}
unsigned int SWIFT::get_m(const HDistribution& distr,ffloat tau){
    unsigned int i=0;
    while(distr.int_error(i++,tau)>TRUNCATION_PRECISION);
    return i;
}
std::complex<ffloat> SWIFT::H(ffloat y, int j){
    ffloat u_j=M_PI*(2*static_cast<ffloat>(j)+1)/(2.*static_cast<ffloat>(J))*static_cast<ffloat>(exp2_m);
    //std::cout<<"u_j"<<u_j<<'\n';
    return 1i*std::exp(1i*u_j*y)*(1/u_j-std::exp(y)/(-1i+u_j));
}
void SWIFT::get_FFT_coeffs(){
    int i;
    fftw_complex * payoff_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *2*J);
    fftw_complex * payoff_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *2*J);
    ffloat * density_in=(ffloat*) fftw_malloc(sizeof(ffloat)*2*J);
    fftw_complex* density_out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(J+1));
    if(!payoff_in||!payoff_out||!density_in||!density_out) std::runtime_error("Wrong FFT params");
    fftw_plan payoff_plan = fftw_plan_dft_1d(2*J, payoff_in, payoff_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan density_plan = fftw_plan_dft_r2c_1d(2*J, density_in, density_out, FFTW_ESTIMATE);
    std::vector<std::complex<ffloat>>& density_coeff_ref=*density_coeffs; //FIXME
    for (i = 0; i < J; ++i)
    {
        std::complex<ffloat> current=H(upper,i)-H(std::max(lower,0.),i);
        payoff_in[i][0] = current.real();
        payoff_in[i][1] = current.imag();
    }
    for (i = J; i < 2*J; ++i)
    {
        payoff_in[i][0] = 0.;
        payoff_in[i][1] = 0.;
    }
    fftw_execute(payoff_plan);
    std::cout<<"First plan executed\n";
    for(i=0;i<2*J;i++) density_in[i]=0.;
    for(i=k_1;i<=k_2;i++){
        //std::cout<<"i mod 2J: "<<i%(2*J)<<"\ti: "<<i<<"\t2J: "<<2*J<<'\n';
        unsigned int i_m_2J=i<0 ? 2*J+i :i;
        std::complex<ffloat> current(payoff_out[i_m_2J][0],payoff_out[i_m_2J][1]);
        //std::cout <<"current: "<<current<<'\n';
        density_in[i_m_2J]=(std::exp(1i*static_cast<ffloat>(i)*M_PI/static_cast<ffloat>(2*J))*current).real()*sqrt_exp2_m/static_cast<ffloat>(J);
    }
    fftw_execute(density_plan);
    std::cout<<"Second plan executed\n";
    for(i=0;i<J/2-1;i++){
        std::complex<ffloat> current(density_out[2*i+1][0],density_out[2*i+1][1]);
        std::cout <<"current: "<<current<<'\n';
        density_coeff_ref[i]=current;
    }
    for(i=J/2-1;i<J;i++){
        std::complex<ffloat> current(density_out[2*J-(2*i+1)][0],density_out[2*J-(2*i+1)][1]);
        std::cout <<"current: "<<current<<'\n';
        density_coeff_ref[i]=current;
    }
    fftw_destroy_plan(payoff_plan);
    fftw_destroy_plan(density_plan);
}
SWIFT::~SWIFT(){
    delete density_coeffs;
}

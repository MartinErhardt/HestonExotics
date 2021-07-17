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

SWIFT::SWIFT(HDistribution& init_distr, const swift_parameters& init_params) : distr(init_distr), my_params(swift_parameters(init_params)), density_coeffs(*(new std::vector<std::complex<ffloat>>(init_params.J))){
    get_FFT_coeffs();
}
std::unique_ptr<swift_parameters> SWIFT::get_parameters(const HDistribution& distr,const ffloat stock_price, const options_chain& opts){
    unsigned int m=0;
    ffloat tau=static_cast<ffloat>(opts.days_to_expiry)/trading_days;
    while(distr.int_error(m++,tau)>TRUNCATION_PRECISION);
    std::cout<<"min_strike: "<<opts.min_strike<<"\tmax_strike: "<<opts.max_strike<<'\n';
    ffloat max=yearly_risk_free*tau+std::log(stock_price/opts.min_strike);
    ffloat min=yearly_risk_free*tau+std::log(stock_price/opts.max_strike);
    ffloat c=std::abs(distr.first_order_moment(tau))+10.*std::sqrt(std::fabs(distr.second_order_moment(tau))+std::sqrt(std::abs(distr.fourth_order_moment(tau))));
    unsigned int exp2_m=std::exp2(m);
    ffloat sqrt_exp2_m=std::sqrt(static_cast<ffloat>(exp2_m));
    ffloat lower=min-c;
    ffloat upper=max+c;
    int k_1=ceil(exp2_m*lower);
    int k_2=floor(exp2_m*upper);
    ffloat iota_density=ceil(std::log2(M_PI*std::abs(k_1-k_2)))-1;
    unsigned int J=std::exp2(iota_density-1);
    std::cout<<"m: "<<m<<"\texp2_m: "<<exp2_m<<"\tc:"<<c<<"\tLower integral bound: "<<k_1<<"\tUpper integral bound: "<<k_2<<"\tJ: "<<J<<'\n';
    return std::unique_ptr<swift_parameters>(new swift_parameters({m,exp2_m,sqrt_exp2_m,lower,upper,k_1,k_2,J}));
}
void SWIFT::get_FFT_coeffs(){
    auto const& [m, exp2_m,sqrt_exp2_m, lower,upper, k_1,k_2,J] = my_params;
    int i;
    unsigned int j;
    auto H = [J,exp2_m](ffloat y, ffloat exp_y, int j)
    {
            ffloat u_j=M_PI*(2*static_cast<ffloat>(j)+1)/(2.*static_cast<ffloat>(J))*static_cast<ffloat>(exp2_m);
            //std::cout<<"u_j"<<u_j<<'\n';
            return 1i*std::exp(1i*u_j*y)*(1/u_j-exp_y/(-1i+u_j));
    };
    ffloat exp_upper=std::exp(upper);
    ffloat exp_lower=std::exp(std::max(lower,0.));
    fftw_complex * payoff_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *2*J);
    fftw_complex * payoff_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *2*J);
    ffloat * density_in=(ffloat*) fftw_malloc(sizeof(ffloat)*4*J);
    fftw_complex* density_out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(2*J+1));
    if(!payoff_in||!payoff_out||!density_in||!density_out) std::runtime_error("Wrong FFT params");
    fftw_plan payoff_plan = fftw_plan_dft_1d(2*J, payoff_in, payoff_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan density_plan = fftw_plan_dft_r2c_1d(4*J, density_in, density_out, FFTW_ESTIMATE);
    for (j = 0; j < J; j++)
    {
        std::complex<ffloat> current=H(upper,exp_upper,j)-H(std::max(lower,0.),exp_lower,j);
        payoff_in[j][0] = current.real();
        payoff_in[j][1] = current.imag();
    }
    for (j = J; j < 2*J; j++)
    {
        payoff_in[j][0] = 0.;
        payoff_in[j][1] = 0.;
    }
    fftw_execute(payoff_plan);
    std::cout<<"First plan executed\n";
    for(j=0;j<4*J;j++) density_in[j]=0.;
    for(i=k_1;i<=k_2;i++){
        //std::cout<<"i mod 2J: "<<i%(2*J)<<"\ti: "<<i<<"\t2J: "<<2*J<<'\n';
        unsigned int i_mod_2J=(2*J+i)&(2*J-1); // This is possible, because J is a power of two.
        unsigned int i_mod_4J=(4*J+i)&(4*J-1);
        //std::complex<ffloat> current(payoff_out[i_m_2J][0],payoff_out[i_m_2J][1]);
        //std::cout <<"current: "<<current<<'\n';
        density_in[i_mod_4J]=(std::exp(1i*static_cast<ffloat>(i)*M_PI/static_cast<ffloat>(2*J))*(payoff_out[i_mod_2J][0]+1i*payoff_out[i_mod_2J][1])).real()*sqrt_exp2_m/static_cast<ffloat>(J);
    }
    fftw_execute(density_plan);
    std::cout<<"Second plan executed\n";
    for(j=0;j<J;j++){
        //std::complex<ffloat> current(density_out[2*i+1][0],density_out[2*i+1][1]);
        std::cout <<"current: "<<density_out[2*j+1][0]-1i*density_out[2*j+1][1]<<'\n';
        density_coeffs[j]=density_out[2*j+1][0]-1i*density_out[2*j+1][1];
    }
    fftw_destroy_plan(payoff_plan);
    fftw_destroy_plan(density_plan);
}
SWIFT::~SWIFT(){
    delete &density_coeffs;
}

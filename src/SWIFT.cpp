/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"SWIFT.h"
#include"my_math.h"
#include <iostream>
#include<fftw3.h>
#include<complex>

using namespace std::complex_literals;

#ifdef FASTMATH
#define TRUNCATION_PRECISION 1e-7
#else
#define TRUNCATION_PRECISION 1e-10
#endif

ffloat SwiftParameters::u(const unsigned int i) const{
    return M_PI*(2*static_cast<ffloat>(i)+1)/(2.*static_cast<ffloat>(J))*static_cast<ffloat>(exp2_m);
}
SwiftParameters::SwiftParameters(const HDistribution& distr,const ffloat stock_price, const options_chain& opts){
    m=0;
    while(distr.int_error(m++)>TRUNCATION_PRECISION); 
    ffloat max=distr.risk_free*distr.tau+std::log(stock_price/opts.min_strike);
    ffloat min=distr.risk_free*distr.tau+std::log(stock_price/opts.max_strike);
    ffloat c=std::abs(distr.first_order_moment())+10.*std::sqrt(std::fabs(distr.second_order_moment())+std::sqrt(std::abs(distr.fourth_order_moment())));
    exp2_m=std::exp2(m);
    sqrt_exp2_m=std::sqrt(static_cast<ffloat>(exp2_m));
    lower=min-c;
    upper=max+c;
    k_1=ceil(exp2_m*lower);
    k_2=floor(exp2_m*upper);
    ffloat iota_density=ceil(std::log2(M_PI*std::abs(k_1-k_2)))-1;
    J=std::exp2(iota_density-1);
}

SWIFT::SWIFT(const swift_parameters& init_params) : my_params(swift_parameters(init_params)), density_coeffs(init_params.J) {
    auto const& [m, exp2_m,sqrt_exp2_m, lower,upper, k_1,k_2,J] = my_params;
    int i;
    unsigned int j;
    auto H = [this](ffloat y, ffloat exp_y, int j){
            return -1i*std::exp(-1i*this->my_params.u(j)*y)*(1/this->my_params.u(j)-exp_y/(1i+this->my_params.u(j)));
    };
    ffloat exp_upper = std::exp(upper);
    ffloat exp_lower = std::exp(std::max(lower,0.));
    fftw_complex * payoff_in =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *2*J);
    fftw_complex * payoff_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *2*J);
    ffloat * density_in =       (ffloat*)       fftw_malloc(sizeof(ffloat)*4*J);
    fftw_complex* density_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(2*J+1));
    if(!payoff_in||!payoff_out||!density_in||!density_out) std::runtime_error("Wrong FFT params!");
    fftw_plan payoff_plan =  fftw_plan_dft_1d(2*J, payoff_in, payoff_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan density_plan = fftw_plan_dft_r2c_1d(4*J, density_in, density_out, FFTW_ESTIMATE);
    //std::cout.precision(dbl::max_digits10);
    for (j = 0; j < J; j++){
        std::complex<ffloat> current=H(upper,exp_upper,j)-H(std::max(lower,0.),exp_lower,j);
        payoff_in[j][0] = current.real();
        payoff_in[j][1] = current.imag();
    }
    for (j = J; j < 2*J; j++){
        payoff_in[j][0] = 0.;
        payoff_in[j][1] = 0.;
    }
    fftw_execute(payoff_plan);
    for(j=0;j<4*J;j++) density_in[j]=0.;
    for(i=k_1;i<=k_2;i++){
        unsigned int i_mod_2J=(2*J+i)&(2*J-1); // This is possible, because J is a power of two.
        unsigned int i_mod_4J=(4*J+i)&(4*J-1);
        std::complex<ffloat> current(payoff_out[i_mod_2J][0],payoff_out[i_mod_2J][1]);
        density_in[i_mod_4J]=(std::exp(1i*static_cast<ffloat>(i)*M_PI/static_cast<ffloat>(2*J))*(payoff_out[i_mod_2J][0]+1i*payoff_out[i_mod_2J][1])).real()*sqrt_exp2_m/static_cast<ffloat>(J);
    }
    fftw_execute(density_plan);
    for(j=0;j<J;j++) density_coeffs[j]=(density_out[2*j+1][0]-1i*density_out[2*j+1][1])*sqrt_exp2_m/static_cast<ffloat>(J); //Why the sqrt_exp2_m/static_cast<ffloat>(J) factor? Also note the minus: We conjugate because we do a backwards transformation, but FFTW does a forward transformation so we read from the back using conjugate symmetry
    fftw_destroy_plan(payoff_plan); // does a fftw_free on input and output
    fftw_destroy_plan(density_plan);
    fftw_free(density_in);          // ... so I thought until I ran address-sanitizer
    fftw_free(density_out);
    fftw_free(payoff_in);
    fftw_free(payoff_out);
}
SWIFT::cache_entry * SWIFT::get_precached(const HDistribution& distr,const ffloat S, const options_chain& opts){
    cache_entry* precached=nullptr; 
    const auto is_opts= [&opts] (const cache_entry& e) {return &e.to_price==&opts;};
    if(!results_cache.size()||!is_opts(*(precached=&*find_if(results_cache.begin(), results_cache.end(), is_opts))))
        precached=&results_cache.emplace_back(distr,*this,opts,S);
    return precached;
}
void SWIFT::price_opts(const HDistribution& distr,const ffloat S, const options_chain& opts,ffloat** out,ffloat* end){
    unsigned int i;
    cache_entry* precached=get_precached(distr,S,opts);
    for(i=0;i<opts.options.size()&&(*out)<end;i++) *((*out)++)=precached->results(0,i).real();
    assert(i>=opts.options.size());
}
void SWIFT::price_opts_grad(const HDistribution& distr,const ffloat S, const options_chain& opts, ffloat** out,ffloat* end){
    unsigned int i;
    cache_entry* precached=get_precached(distr,S,opts);
    for(i=0;i<opts.options.size()&&(*out)<end;i++) for(unsigned int j=1;j<6;j++) *((*out)++)=precached->results(j,i).real();
    assert(i>=opts.options.size());
}
void SWIFT::flush_cache(){
    results_cache.clear();
}
SWIFT::CacheEntry::CacheEntry(const HDistribution& distr, const SWIFT& swift_obj, const options_chain& to_price_init,const ffloat stock_price):to_price(to_price_init){
    unsigned int swift_J=swift_obj.my_params.J;
    ffloat discount=std::exp(-distr.risk_free*distr.tau);
    MatC pricing_matrix(6,swift_J);
    MatC to_price_matrix(swift_J,to_price.options.size());
    for(unsigned int i=0;i<swift_J;i++){
        std::vector<std::complex<ffloat>> chf_chf_grad_val=distr.chf_chf_grad(swift_obj.my_params.u(i));
        for(int j=0;j<6;j++) pricing_matrix(j,i)=chf_chf_grad_val[j]*swift_obj.density_coeffs[i];
    }
    unsigned int j=0;
    for(const option& opt: to_price.options){
        ffloat x=distr.risk_free*distr.tau+std::log(stock_price/opt.strike);
        for(unsigned int i=0;i<swift_J;i++) to_price_matrix(i,j)=discount*opt.strike*std::exp(-1i*swift_obj.my_params.u(i)*x);
        j++;
    }
    results=pricing_matrix*to_price_matrix;
}
std::ostream& operator<<(std::ostream& out, SwiftParameters const& sp)
{
    return out <<"m: "<< sp.m << ", " << "exp2_m: "<< sp.exp2_m << ", " << "sqrt_exp2_m: "<< sp.sqrt_exp2_m << ", "
               <<"k1: "<< sp.k_1 << ", " << "lower: "<< sp.lower << ", " << "k2: "<< sp.k_2 << ", " << "upper: "<< sp.upper << ", "
               << "J: "<< sp.J << std::endl;
}

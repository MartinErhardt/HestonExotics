/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"SWIFT.h"
#include <iostream>
#include<fftw3.h>
#include<complex>

using namespace std::complex_literals;
//typedef Matrix<std::complex<ffloat>, 6, Eigen::Dynamic> MatrixXdcf6;
//typedef Matrix<std::complex<ffloat>, Eigen::Dynamic, Eigen::Dynamic> MatrixXdcf;
#ifdef FASTMATH
#define TRUNCATION_PRECISION 1e-7
#else
#define TRUNCATION_PRECISION 1e-10
#endif

ffloat SwiftParameters::u(const unsigned int i) const{
    return M_PI*(2*static_cast<ffloat>(i)+1)/(2.*static_cast<ffloat>(J))*static_cast<ffloat>(exp2_m);
}
std::unique_ptr<swift_parameters> SWIFT::get_parameters(const HDistribution& distr,const ffloat stock_price, const options_chain& opts){
    //std::cout<<"tau: "<<distr.tau<<"\tmin_strike: "<<opts.min_strike<<"\tmax_strike: "<<opts.max_strike<<'\n';
    unsigned int m=0;
    while(distr.int_error(m++)>TRUNCATION_PRECISION); 
    //std::cout<<"min_strike: "<<opts.min_strike<<"\tmax_strike: "<<opts.max_strike<<'\n';
    ffloat max=distr.risk_free*distr.tau+std::log(stock_price/opts.min_strike);
    ffloat min=distr.risk_free*distr.tau+std::log(stock_price/opts.max_strike);
    ffloat c=std::abs(distr.first_order_moment())+10.*std::sqrt(std::fabs(distr.second_order_moment())+std::sqrt(std::abs(distr.fourth_order_moment())));
    unsigned int exp2_m=std::exp2(m);
    ffloat sqrt_exp2_m=std::sqrt(static_cast<ffloat>(exp2_m));
    ffloat lower=min-c;
    ffloat upper=max+c;
    int k_1=ceil(exp2_m*lower);
    int k_2=floor(exp2_m*upper);
    ffloat iota_density=ceil(std::log2(M_PI*std::abs(k_1-k_2)))-1;
    unsigned int J=std::exp2(iota_density-1);
    //std::cout<<"m: "<<m<<"\texp2_m: "<<exp2_m<<"\tc:"<<c<<"\tLower integral bound: "<<k_1<<"\tUpper integral bound: "<<k_2<<"\tJ: "<<J<<'\n';
    return std::unique_ptr<swift_parameters>(new swift_parameters({m,exp2_m,sqrt_exp2_m,lower,upper,k_1,k_2,J}));
}

SWIFT::SWIFT(const swift_parameters& init_params) : my_params(swift_parameters(init_params)), density_coeffs(*(new std::vector<std::complex<ffloat>>(init_params.J))),results_cache(*(new std::list<cache_entry>())) {
    get_FFT_coeffs();
}
void SWIFT::get_FFT_coeffs(){
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
        //std::cout<<"j: "<<j<<"\tval: "<<current<<'\n';
    }
    for (j = J; j < 2*J; j++){
        payoff_in[j][0] = 0.;
        payoff_in[j][1] = 0.;
    }
    fftw_execute(payoff_plan);
    //std::cout<<"First plan executed\n";
    for(j=0;j<4*J;j++) density_in[j]=0.;
    for(i=k_1;i<=k_2;i++){
        //std::cout<<"i mod 2J: "<<i%(2*J)<<"\ti: "<<i<<"\t2J: "<<2*J<<'\n';
        unsigned int i_mod_2J=(2*J+i)&(2*J-1); // This is possible, because J is a power of two.
        unsigned int i_mod_4J=(4*J+i)&(4*J-1);
        std::complex<ffloat> current(payoff_out[i_mod_2J][0],payoff_out[i_mod_2J][1]);
        //std::cout <<"i: "<<i<<"\ti_mod_N: "<<i_mod_2J<<"\tidfpart: "<<current<<'\n';
        density_in[i_mod_4J]=(std::exp(1i*static_cast<ffloat>(i)*M_PI/static_cast<ffloat>(2*J))*(payoff_out[i_mod_2J][0]+1i*payoff_out[i_mod_2J][1])).real()*sqrt_exp2_m/static_cast<ffloat>(J);
    }
    fftw_execute(density_plan);
    //std::cout<<"Second plan executed\n";
    for(j=0;j<J;j++){
        //std::complex<ffloat> current(density_out[2*i+1][0],density_out[2*i+1][1]);
        density_coeffs[j]=(density_out[2*j+1][0]-1i*density_out[2*j+1][1])*sqrt_exp2_m/static_cast<ffloat>(J); //Why the sqrt_exp2_m/static_cast<ffloat>(J) factor?
        //std::cout <<"u_"<<j<<": "<<density_coeffs[j]<<'\n';
    }
    fftw_destroy_plan(payoff_plan);
    fftw_destroy_plan(density_plan);
}
SWIFT::cache_entry * SWIFT::get_precached(const HDistribution& distr,const ffloat S, const options_chain& opts){
    cache_entry* precached=nullptr; 
    const auto is_opts= [&opts] (const cache_entry& e) {return &e.to_price==&opts;};
    if(!results_cache.size()||!is_opts(*(precached=&*find_if(results_cache.begin(), results_cache.end(), is_opts)))) //TODO test if this &* is local in scope
        precached=&results_cache.emplace_back(distr,*this,opts,S);
    //else std::cout<<"Used precached values\n";
    //std::cout<<"returning precached or newly cached results\t# rows: "<<precached->results.rows()<<"\t# cols: "<<precached->results.cols()<<'\n';
    return precached;
}
void SWIFT::price_opts(const HDistribution& distr,const ffloat S, const options_chain& opts,ffloat** out,ffloat* end){
    unsigned int i;
    cache_entry* precached=get_precached(distr,S,opts);
    for(i=0;i<opts.options->size()&&(*out)<end;i++){
        *((*out)++)=precached->results(0,i).real();
        //std::cout<<"S: "<<S<<"\tstrike: "<<(*opts.options)[i].strike<<"\tvolume: "<<(*opts.options)[i].volume<<"\tp: "<<*((*out)-1)<<"\treal p: "<<(*opts.options)[i].price<<"\tdiff: "<<*((*out)-1)-(*opts.options)[i].price<<"\texpiry: "<<opts.time_to_expiry<<'\n';
    }
    if(i<opts.options->size()) throw std::runtime_error((std::string("Pricing buffer too small i: ")+ std::to_string(i)+std::string("\t# strikes: ")+std::to_string(opts.options->size())).c_str());
}
void SWIFT::price_opts_grad(const HDistribution& distr,const ffloat S, const options_chain& opts, ffloat** out,ffloat* end){
    unsigned int i;
    cache_entry* precached=get_precached(distr,S,opts);
    for(i=0;i<opts.options->size()&&(*out)<end;i++) for(unsigned int j=1;j<6;j++){*((*out)++)=precached->results(j,i).real();
        //std::cout<<"grad_"<<j<<": "<<*((*out)-1)<<'\n';
    }
    if(i<opts.options->size()) throw std::runtime_error((std::string("Gradient buffer too small i: ")+std::to_string(i)+std::string("\t# strikes: ")+std::to_string(opts.options->size())).c_str());
}
void SWIFT::flush_cache(){
    results_cache.clear();
}
SWIFT::CacheEntry::CacheEntry(const HDistribution& distr, const SWIFT& swift_obj, const options_chain& to_price_init,const ffloat stock_price):to_price(to_price_init){
    unsigned int swift_J=swift_obj.my_params.J;
    //ffloat e_m=static_cast<ffloat>(swift_obj.my_params.exp2_m);
    ffloat discount=std::exp(-distr.risk_free*distr.tau);
    //std::cout<<"J: "<<swift_J<<"\t# obs: "<<to_price.options->size();
    MatC pricing_matrix(6,swift_J);
    MatC to_price_matrix(swift_J,to_price.options->size());
    //std::cout<<"J: "<<swift_J<<"\tn opts: "<<to_price.options.size()<<'\n';
    //std::cout<<"I allocate "<<sizeof(to_price_matrix)+sizeof(pricing_matrix)<<" bytes on the stack\n";
    for(unsigned int i=0;i<swift_J;i++){
        std::vector<std::complex<ffloat>> chf_chf_grad_val=distr.chf_chf_grad(swift_obj.my_params.u(i));
        
        //std::cout<<"j"<<i<<"u: "<<swift_obj.my_params.u(i)
        //<<"\tchf: "<<chf_chf_grad_val[0]<<"\tv0: "<<distr.p.v_0<<"\tv_m: "<<distr.p.v_m<<"\trho: "<<distr.p.rho<<"\tkappa: "<<distr.p.kappa<<"\tsigma: "<<distr.p.kappa<<'\n';
        
        //std::cout<<"i: "<<i<<"\ta_i0: "<<chf_chf_grad_val[i]<<'\n';
        for(int j=0;j<6;j++){ pricing_matrix(j,i)=chf_chf_grad_val[j]*swift_obj.density_coeffs[i];
            //std::cout<<"i: "<<i<<"\tj: "<<j<<"\ttau: "<<distr.tau<<"\tu_i: "<<swift_obj.my_params.u(i)<<"\tval: "<<chf_chf_grad_val[j]<<'\n';
            
        }
    }
    unsigned int j=0;
    //std::cout<<"Filled pricing matrix\n";
    for(const option& opt: *to_price.options){
        ffloat x=distr.risk_free*distr.tau+std::log(stock_price/opt.strike);
        for(unsigned int i=0;i<swift_J;i++){ to_price_matrix(i,j)=discount*opt.strike*std::exp(-1i*swift_obj.my_params.u(i)*x);
        //std::cout<<"j: "<<j<<"\ti: "<<i<<"\tb_ji: "<<discount*opt.strike*std::exp(-1i*swift_obj.my_params.u(i)*x)<<'\n';
            
        }
        j++;
    }
    //std::cout<<"Filled to price matrix\n";
    results=pricing_matrix*to_price_matrix;
}

SWIFT::~SWIFT(){
    //std::cout<<"delete";
    delete &density_coeffs;
    delete &results_cache;
}

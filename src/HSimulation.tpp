/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HSimulation.h"
#include<numeric>
#include"RNG.h"
#include <omp.h>
using namespace HSimulation;

template<class Scheme>
std::vector<ffloat> HSimulation::price(const HParams& p, const ffloat S, 
                                        const std::list<options_chain>& all_chains,
                                        unsigned int n_simulations, unsigned int n_opts,unsigned int steps){
    // assert, that expiration dates are in increasing order
#ifndef NDEBUG
    typedef struct {bool exp_date_increasing; ffloat prev_exp_date;} assert_increasing_s;
#endif
    assert(std::accumulate(std::next(all_chains.begin()), all_chains.end(), assert_increasing_s({true,all_chains.begin()->time_to_expiry}), 
                                [](assert_increasing_s cum, auto const& it) {
                                    return assert_increasing_s({cum.exp_date_increasing&&cum.prev_exp_date<it.time_to_expiry,it.time_to_expiry});
                                }).exp_date_increasing);
    std::vector<ffloat> prices=*(new std::vector<ffloat>(n_opts));
#pragma omp parallel
{
    std::vector<ffloat> local_prices=*(new std::vector<ffloat>(n_opts));
    const SDE_state<2> initial_state={{S,p.v_0},{S,p.v_0},.0,.0};
    const unsigned int local_sims=n_simulations/omp_get_num_threads();
    RNG local_rng(RAND_BUF_SIZE,(1<<omp_get_thread_num()));
    Scheme heston_sde(p,all_chains.begin()->time_to_expiry/steps,&local_rng,initial_state);
    for(unsigned int i=0; i<local_sims;i++){
        unsigned int opts_priced=0;
        std::list<options_chain>::const_iterator earliest_unpriced=all_chains.begin();
        heston_sde.update_earliest(earliest_unpriced->time_to_expiry,steps);
        heston_sde.reset();
        for(++(heston_sde=initial_state);opts_priced<n_opts;++heston_sde){
            while(heston_sde>=earliest_unpriced->time_to_expiry&&opts_priced<n_opts){
                unsigned int opt_index=0;
                heston_sde.accumulate_final_value(heston_sde.state);
                for(const auto& opt:earliest_unpriced->options)
                    prices[opts_priced+(opt_index++)]+=heston_sde.final_payoff(opt.strike)/n_simulations;
                if((opts_priced=opts_priced+earliest_unpriced->options.size())<n_opts)
                    heston_sde.update_earliest((++earliest_unpriced)->time_to_expiry,steps);
            }
            heston_sde.accumulate_value(heston_sde.state);
        }
    }
#pragma omp critical
    for(unsigned int i=0;i<n_opts;i++) prices[i]+=local_prices[i]/omp_get_num_threads();
}
    return prices;
}
template<typename accumulate_t,class OptionPolicy> 
SDE<2>& HQEAnderson<accumulate_t,OptionPolicy>::operator++(){
    auto const& [v_0,theta,rho,kappa,eps] = params;
    ffloat delta=step_width(state);
    ffloat gamma_1=.5;
    ffloat gamma_2=.5;
    ffloat discount=std::exp(-kappa*delta);
    ffloat m=theta+(state.cur[V]-theta)*discount;
    ffloat sp2=std::fabs(state.cur[V]*eps*eps*discount/kappa*(1-discount)+theta*eps*eps/(2*kappa)*(1.-discount)*(1.-discount));
    ffloat Psi=sp2/(m*m);
    state.prev[V]=state.cur[V];
    if(Psi<PSI_C){
        ffloat bp2=2/Psi-1+std::sqrt(2/Psi*(2/Psi-1));
        ffloat b=std::sqrt(bp2);
        ffloat a=m/(1+bp2);
        ffloat Z_V=rng->get_grand();
        state.cur[V]=a*(b+Z_V)*(b+Z_V);
    } else{
        ffloat p=(Psi-1)/(Psi+1);
        ffloat beta=2/(m*(Psi+1));
        ffloat U_V=rng->get_urand();
        state.cur[V]=p<U_V? std::log((1-p)/(1-U_V))/beta:0.;
    }
    ffloat K_0=-rho*kappa*theta/eps*delta;
    ffloat K_1=gamma_1*delta*(kappa*rho/eps-.5)-rho/eps;
    ffloat K_2=gamma_2*delta*(kappa*rho/eps-.5)+rho/eps;
    ffloat K_3=gamma_1*delta*(1-rho*rho);
    ffloat K_4=gamma_2*delta*(1-rho*rho);
    log_X=log_X+K_0+K_1*state.prev[V]+K_2*state.cur[V]+std::sqrt(K_3*state.prev[V]+K_4*state.cur[V])*rng->get_grand();
    state.prev[X]=state.cur[X];
    state.cur[X]=std::exp(log_X);
    state.prev_time=state.cur_time;
    state.cur_time+=delta;
    return *this;
}
template<typename accumulate_t,class OptionPolicy>
HQEAnderson<accumulate_t,OptionPolicy>& HQEAnderson<accumulate_t,OptionPolicy>::operator=(const SDE_state<2> new_state) {
    state=new_state;
    log_X=std::log(state.cur[X]);
    state.cur_time=new_state.cur_time;
    state.prev_time=new_state.prev_time;
    return *this;
}

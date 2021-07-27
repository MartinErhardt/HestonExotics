/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
//#include <limits>
#include"AsianContract.h"
#include"HSimulation.h"
using namespace HSimulation;
//using HQEAnderson<ffloat,NonAdaptive>::positions;
std::vector<ffloat>* ArithmeticAOptionTool::price_trapezoidal(const HParams& p, const ffloat s, 
                                                              const std::vector<options_chain>& all_chains,
                                                              unsigned int n_opts,unsigned int min_steps){
    const SDE_state<2> initial_state={{s,p.v_0},{0.0,0.0},0.0,0.0};
    HQEAnderson<ffloat,NonAdaptive> heston_sde(p,all_chains.begin()->time_to_expiry/min_steps,my_rng,initial_state);
    std::vector<ffloat>& prices=*(new std::vector<ffloat>(n_opts));
    //TODO assert all_chains ordered by expiry
    //TODO assert n_opts # of all options
    for(unsigned int i=0; i<n_simulations;i++){
        ffloat underlying_integral=0.;
        unsigned int opts_priced=0;
        std::vector<options_chain>::const_iterator earliest_unpriced=all_chains.begin();
        heston_sde.init_step_size=earliest_unpriced->time_to_expiry/min_steps;
        for(++(heston_sde=initial_state);opts_priced<n_opts;++heston_sde){
            while(heston_sde>=earliest_unpriced->time_to_expiry&&opts_priced<n_opts){
                unsigned int opt_index=0;
                ffloat step_interpolation=(heston_sde->cur[X]-heston_sde->prev[X])*(earliest_unpriced->time_to_expiry-heston_sde->prev_time)
                                                    /heston_sde.init_step_size;
                ffloat underlying_mean=(underlying_integral+step_interpolation)/earliest_unpriced->time_to_expiry;
                for(const auto& opt:*earliest_unpriced->options)
                    prices[opts_priced+opt_index++]+=std::max(underlying_mean-opt.strike,0.)/n_simulations;
                if(opts_priced+=earliest_unpriced->options->size()<n_opts)
                    heston_sde.init_step_size=(++earliest_unpriced)->time_to_expiry/min_steps;
            }
            underlying_integral+=heston_sde.init_step_size*.5*(heston_sde->cur[X]+heston_sde->prev[X]); //TODO enum
        }
    }
    return &prices;
}

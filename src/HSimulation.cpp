/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HSimulation.h"
using namespace HSimulation;

template<typename accumulate_t,class Scheme>
std::vector<ffloat>* PricingTool<accumulate_t,Scheme>::price(const HParams& p, const ffloat S, 
                                        const std::list<options_chain>& all_chains,
                                        unsigned int n_simulations, unsigned int n_opts,unsigned int steps){
    const SDE_state<2> initial_state={{S,p.v_0},{S,p.v_0},.0,.0};
    HQEAnderson<accumulate_t,Scheme> heston_sde(p,all_chains.begin()->time_to_expiry/steps,my_rng,initial_state);
    std::vector<ffloat>& prices=*(new std::vector<ffloat>(n_opts));
    //TODO assert all_chains ordered by expiry
    //TODO assert n_opts # of all options
    //ffloat latest_unpriced=all_chains.end()->time_to_expiry;
    for(unsigned int i=0; i<n_simulations;i++){
        unsigned int opts_priced=0;
        std::list<options_chain>::const_iterator earliest_unpriced=all_chains.begin();
        heston_sde.update_earliest(earliest_unpriced->time_to_expiry,steps);
        heston_sde.reset();
        //for(const auto& opt:*earliest_unpriced->options)
        for(++(heston_sde=initial_state);opts_priced<n_opts;++heston_sde){
            while(heston_sde>=earliest_unpriced->time_to_expiry&&opts_priced<n_opts){
                unsigned int opt_index=0;
                heston_sde.accumulate_final_value(heston_sde.state);
                for(const auto& opt:*earliest_unpriced->options)
                    prices[opts_priced+(opt_index++)]+=heston_sde.final_payoff(opt.strike)/n_simulations;
                if((opts_priced=opts_priced+earliest_unpriced->options->size())<n_opts)
                    heston_sde.update_earliest((++earliest_unpriced)->time_to_expiry,steps);
            }
            heston_sde.accumulate_value(heston_sde.state);
        }
    }
    return &prices;
}

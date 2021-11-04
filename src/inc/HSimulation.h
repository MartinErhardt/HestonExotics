/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include<memory>
#include <list>
#include"SDE.h"
#include"Types.h"
#include"HDistribution.h"
#include"RNG.h"
#define RAND_BUF_SIZE (1<<20)//16 //128/8
#define PSI_C 1.5 ///< point at which the process is so small, that it makes sense to switch methods

namespace HSimulation{
    /**enum containing positions of price process and volatility process in the two-dimensional SDE*/
    enum positions{X=0,V=1};
    /**
     * @brief template to simulate the Heston model SDE using the Quadratic Exponential (QE) method introduced by Anderson.
     * 
     * This uses an approach called Policy-based design. Scheme is another class that contains a function called step_width. Using this function we can avoid implementing any function step_width and customize step_width as if it was a virtual function. Furthermore this leads to a much better performance.
     */
    template<typename accumulate_t,class OptionPolicy> class HQEAnderson : OptionPolicy,SDE<2>{
        ffloat log_X;
        const HParams params;    
    public:
        using OptionPolicy::step_width;
        using OptionPolicy::accumulate_value;
        using OptionPolicy::accumulate_final_value;
        using OptionPolicy::reset;
        using OptionPolicy::update_earliest;
        using OptionPolicy::final_payoff;
        using SDE::operator->;
        using SDE::operator*;
        using SDE::operator>=;
        using SDE::operator++;
        using SDE::state;
        HQEAnderson(const HParams& params_init,accumulate_t option_policy_init, RNG*rng, const SDE_state<2> init_cond) : OptionPolicy(option_policy_init),SDE<2>(init_cond,rng),log_X(std::log(init_cond.cur[0])),params(params_init){}
        /**
         * @brief implements a step of the Quadratic Exponential (QE) method to simulate the Heston SDE, first described by Anderson
         * 
         * Note that this has to be implemented inline.
         * 
         * @return updated this instance
         */
        SDE<2>& operator++();
        HQEAnderson<accumulate_t,OptionPolicy>& operator=(const SDE_state<2> new_state);
    };
    /**
    * @brief prices all options in all of the options chains in all_chains using the trapezoidal rule
    * 
    * This could use a control variate approach, that utilizes the price of a geometric asian option as such a control variate.
    * 
    * @param p Heston model parameters
    * @param S stock price
    * @param all_chains vector of options chains
    * @param n_observations total number of options to be priced
    * @param min_steps minimum number of step for every any option to be computed
    * @return pointer to std::vector of option prices ordered in the the same ways as the options of all_chains and chain.options respectively
    * @throws abort if NDEBUG is enabled and all_chains are not in increasing order by expiry date
    */
    template<class Scheme>
    std::vector<ffloat>* price(const HParams& p, const ffloat S,const std::list<options_chain>& all_chains,
                                            unsigned int n_simulations, unsigned int n_opts,unsigned int min_steps);
}

#include "../HSimulation.tpp"

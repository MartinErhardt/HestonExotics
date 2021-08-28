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
    template<typename accumulate_t,class Scheme> class HQEAnderson : Scheme,SDE<2>{
        ffloat log_X;
        const HParams params;    
    public:
        using Scheme::step_width;
        using Scheme::accumulate_value;
        using Scheme::accumulate_final_value;
        using Scheme::update_earliest;
        using Scheme::final_payoff;
        using SDE::operator->;
        using SDE::operator*;
        using SDE::operator>=;
        using SDE::operator++;
        using SDE::state;
        HQEAnderson(const HParams& params_init,accumulate_t scheme_params_init, RNG*rng, const SDE_state<2> init_cond) : Scheme(scheme_params_init),SDE<2>(init_cond,rng),log_X(std::log(init_cond.cur[0])),params(params_init){}
        /**
         * @brief implements a step of the Quadratic Exponential (QE) method to simulate the Heston SDE, first described by Anderson
         * 
         * Note that this has to be implemented inline.
         * 
         * @return updated this instance
         */
        SDE<2>& operator++(){
            auto const& [v_0,theta,rho,kappa,eps] = params;
            ffloat delta=step_width(state);
            ffloat gamma_1=.5;
            ffloat gamma_2=.5;
            ffloat discount=std::exp(-kappa*delta);
            ffloat m=theta+(state.cur[V]-theta)*discount;
            ffloat sp2=state.cur[V]*eps*eps*discount/kappa*(1-discount)+theta*eps*eps/(2*kappa)*(1.-discount)*(1.-discount);
            ffloat Psi=sp2/(m*m);
            ffloat p=(Psi-1)/(Psi+1);
            ffloat beta=2/(m*(Psi+1));
            ffloat bp2=2/Psi-1+std::sqrt(2/Psi*(2/Psi-1));
            ffloat b=std::sqrt(bp2);
            ffloat a=m/(1+bp2);
            state.prev[V]=state.cur[V];
            if(Psi>PSI_C){
                ffloat Z_V=rng->get_grand();
                state.cur[V]=a*(b+Z_V)*(b+Z_V);
            } else{
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
        HQEAnderson<accumulate_t,Scheme>& operator=(const SDE_state<2> new_state) {state=new_state; return *this;}
    };
    /**
     * @brief template to simulate the Heston model SDE using the Quadratic Exponential (QE) method introduced by Anderson.
     * 
     * This uses an approach called Policy-based design. Scheme is another class that contains a function called step_width. Using this function we can avoid implementing any function step_width and customize step_width as if it was a virtual function. Furthermore this leads to a much better performance.
     */
    template<typename accumulate_t, class Scheme>class PricingTool{
        RNG* my_rng;
    public:
        PricingTool(const unsigned int thread_i):my_rng(new RNG(RAND_BUF_SIZE,(1<<thread_i))){};
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
        std::vector<ffloat>* price(const HParams& p, const ffloat S,const std::list<options_chain>& all_chains,
                                                unsigned int n_simulations, unsigned int n_opts,unsigned int min_steps);
        ~PricingTool(){delete my_rng;};
    };
}


/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"
#include"RNG.h"
#include"HDistribution.h"
#define RAND_BUF_SIZE 16 //128/8
/**
 * @brief This class encapsulates a PRNG and a pricing method for arithmetic asian options
 */
class ArithmeticAOptionTool{
    RNG* my_rng;
    const unsigned int n_simulations;
public:
    ArithmeticAOptionTool(const unsigned int thread_n, const unsigned int n_sim_init):my_rng(new RNG(RAND_BUF_SIZE,(1<<thread_n))),n_simulations(n_sim_init){};
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
    std::vector<ffloat>* price_trapezoidal(const HParams& p, const ffloat S,const std::vector<options_chain>& all_chains,unsigned int n_observations,unsigned int min_steps);
    ~ArithmeticAOptionTool(){delete my_rng;};
};

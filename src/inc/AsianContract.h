/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"
#include"RNG.h"
#include"HDistribution.h"
#define RAND_BUF_SIZE 16 //128/8
class ArithmeticAOptionTool{
    RNG* my_rng;
    const unsigned int n_simulations;
public:
    ArithmeticAOptionTool(const unsigned int thread_n, const unsigned int n_sim_init):my_rng(new RNG(RAND_BUF_SIZE,(1<<thread_n))),n_simulations(n_sim_init){};
    std::vector<ffloat>* price_trapezoidal(const HParams& p, const ffloat S,const std::vector<options_chain>& all_chains,unsigned int n_observations,unsigned int min_steps);
    ~ArithmeticAOptionTool(){delete my_rng;};
};

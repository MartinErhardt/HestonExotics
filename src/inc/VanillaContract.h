/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"

#include"HDistribution.h"
#include"HSimulation.h"
#include<iostream>

/**
 * @brief This class specifies the integration scheme to prise a European Call option NonAdaptively
 */
class EuropeanCallNonAdaptive{
    ffloat final_value;
    ffloat earliest_unpriced_expi;
    ffloat init_step_size;
public:
    EuropeanCallNonAdaptive(ffloat fixed_step_size): init_step_size(fixed_step_size){}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
    ffloat step_width(const SDE_state<2>& state){return init_step_size;}
//#pragma GCC diagnostic pop
    void accumulate_value(const SDE_state<2>& state){
//        std::cout<<"t: "<<state.cur_time<<"\tprev t: "<<state.prev_time<<"\tX: "<<state.cur[HSimulation::X]<<"\tV:"<<state.cur[HSimulation::V]<<std::endl;
    }
#pragma GCC diagnostic pop
    void accumulate_final_value(const SDE_state<2>& state){
        final_value=state.prev[HSimulation::X]+(state.cur[HSimulation::X]-state.prev[HSimulation::X])*(earliest_unpriced_expi-state.prev_time)
                                                                /init_step_size;
    }
    void update_earliest(ffloat earliest_expiry,ffloat min_strike){
        init_step_size=earliest_expiry/min_strike;
        earliest_unpriced_expi=earliest_expiry;
    }
    void reset(){}
    ffloat final_payoff(ffloat strike){
        return std::max(0.,final_value-strike);
    }
};
extern HSimulation::HQEAnderson<ffloat,EuropeanCallNonAdaptive> explicit_HQEAnderson_EuropeanCallNonAdaptive;
extern HSimulation::PricingTool<ffloat,EuropeanCallNonAdaptive> explicit_PricingTool_EuropeanCallNonAdaptive;

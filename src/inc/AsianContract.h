/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"

#include"HDistribution.h"
#include"HSimulation.h"
#include<iostream>

/**
 * @brief This class specifies the integration scheme to prise a Arithmetic Asian Call option NonAdaptively
 */
class AAsianCallNonAdaptive{
    ffloat accumulated_value;
    ffloat earliest_unpriced_expi;
    ffloat init_step_size;
public:
    AAsianCallNonAdaptive(ffloat fixed_step_size): init_step_size(fixed_step_size){}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
    ffloat step_width(const SDE_state<2>& state){return init_step_size;}
#pragma GCC diagnostic pop
    void accumulate_value(const SDE_state<2>& state){
        accumulated_value+=init_step_size*.5*(state.cur[HSimulation::X]+state.prev[HSimulation::X]);
        std::cout<<"t: "<<state.cur_time<<"\tprev t: "<<state.prev_time<<"\tX: "<<state.cur[HSimulation::X]<<"\tV:"<<state.cur[HSimulation::V]<<std::endl;
    }
    void accumulate_final_value(const SDE_state<2>& state){
        if(earliest_unpriced_expi>state.cur_time) throw std::runtime_error("not the final value!");
        ffloat step_interpolation=(state.cur[HSimulation::X]-state.prev[HSimulation::X])*(earliest_unpriced_expi-state.prev_time)
                                           /init_step_size;
        accumulated_value=(accumulated_value+step_interpolation)/earliest_unpriced_expi;
    }
    void update_earliest(ffloat earliest_expiry,ffloat min_strike){
        init_step_size=earliest_expiry/min_strike;
        earliest_unpriced_expi=earliest_expiry;
    }
    ffloat final_payoff(ffloat strike){
        return std::max(0.,accumulated_value-strike);
    }
};
extern HSimulation::HQEAnderson<ffloat,AAsianCallNonAdaptive> explicit_HQEAnderson_AAsianCallNonAdaptive;
extern HSimulation::PricingTool<ffloat,AAsianCallNonAdaptive> explicit_PricingTool_AAsianCallNonAdaptive;

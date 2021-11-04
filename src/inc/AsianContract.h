/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"

#include"HDistribution.h"
#include"HSimulation.h"
#include<iostream>
#include"BSM.h"
/**
 * @brief This class specifies the integration scheme to prise a Arithmetic Asian Call option NonAdaptively
 */
class AAsianCallNonAdaptive{
    ffloat accumulated_value;
    ffloat final_value;
    ffloat earliest_unpriced_expi;
    ffloat init_step_size;
public:
    AAsianCallNonAdaptive(ffloat fixed_step_size): init_step_size(fixed_step_size){reset();}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
    ffloat step_width(const SDE_state<2>& state){return init_step_size;}
#pragma GCC diagnostic pop
    void accumulate_value(const SDE_state<2>& state){
        accumulated_value+=init_step_size*.5*(state.cur[HSimulation::X]+state.prev[HSimulation::X]);
        //std::cout<<"t: "<<state.cur_time<<"\tprev t: "<<state.prev_time<<"\tX: "<<state.cur[HSimulation::X]<<"\tprev X: "<<state.prev[HSimulation::X]<<"\tV:"<<state.cur[HSimulation::V]<<"\tprev V: "<<state.prev[HSimulation::V]<<"\tinit_step_size: "<<init_step_size<<"\taccumulated value: "<<accumulated_value<<std::endl;
    }
    void accumulate_final_value(const SDE_state<2>& state){
        if(earliest_unpriced_expi>state.cur_time) throw std::runtime_error("not the final value!");
        ffloat step_interpolation=(state.cur[HSimulation::X]-state.prev[HSimulation::X])*(earliest_unpriced_expi-state.prev_time)
                                                    /init_step_size;
        final_value=(accumulated_value+step_interpolation)/earliest_unpriced_expi;
    }
    void update_earliest(ffloat earliest_expiry,ffloat min_strike){
        init_step_size=earliest_expiry/min_strike;
        earliest_unpriced_expi=earliest_expiry;
    }
    void reset(){
        accumulated_value=0.;
        final_value=0.;
    }
    ffloat final_payoff(ffloat strike){
        ffloat ret_val=std::max(0.,final_value-strike);
        //std::cout<<"strike: "<<strike<<"\tdays to expi: "<<earliest_unpriced_expi*trading_days<<"\tfinal val: "<<final_value<<"\tfinal payoff: "<<ret_val<<std::endl;
        return ret_val;
    }
};
extern HSimulation::HQEAnderson<ffloat,AAsianCallNonAdaptive> explicit_HQEAnderson_AAsianCallNonAdaptive;

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include "RNG.h"
/**
 * @brief template encapsulating the state of a (discretized) SDE of dimension d
 */
template<const unsigned int d> struct SDE_state {
    ffloat cur[d];      //<current SDE values
    ffloat prev[d];     //<previous SDE values (before discretized time step)
    ffloat cur_time;    //<current time
    ffloat prev_time;   //<previous time (before discretized time step)
};
/**
 * @brief This is a extremely thin layer of abstraction simply there to encode iterator propertis and initialize the prng.
 */
template<const unsigned int d> class SDE {
        //using SchemeEvoPolicy::evolution;
    protected:
        RNG * rng;
    public:
        SDE_state<d> state;
        SDE(const SDE_state<d>& init_cond, RNG*init_rng):rng(init_rng),state(init_cond){}
        //SDE & operator++() { evolution(this->state,(rng->buf_cur+=d),step_width(this->state)); return *this;}
        //SDE& operator=(const SDE_state<d> new_state) {state=new_state; return *this;}
        SDE<d>& operator++(int) {SDE& retval = *this; ++(*this); return retval;}
        //bool operator==(SDE other) const {return *this == other;}
        bool operator>=(SDE& other) const {return state.cur_time >= other.state.cur_time;}
        bool operator>=(ffloat time) const {return state.cur_time >= time;}
        //bool operator!=(SDE other) const {return !(*this == other);}
        SDE_state<d>& operator*() {return state;}
        SDE_state<d>* operator->(){return &state;}
        using difference_type = ffloat;
        using value_type = SDE<d>;
        using pointer = const SDE<d> **;
        using reference = const SDE<d> *;
        using iterator_category = std::output_iterator_tag;
};

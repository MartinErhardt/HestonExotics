/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include "RNG.h"
template<const unsigned int d> struct SDE_state {
    double cur_state[d];
    double prev_state[d];
    double time;
};

template<const unsigned int d,typename SchemeEvoParams, class SchemeEvoPolicy> class SDE : private SchemeEvoPolicy{
        SDE_state<d> state;
        RNG * rng;
        using SchemeEvoPolicy::evolution;
        using SchemeEvoPolicy::step_width;
    public:
        SDE(SDE_state<d> * init_cond, RNG * rng,SchemeEvoParams params);
        SDE & operator++() { evolution(this->state,(rng->buf_cur+=d),step_width(this->state)); return *this;}
        SDE operator++(int) {SDE retval = *this; ++(*this); return retval;}
        //bool operator==(SDE other) const {return *this == other;}
        bool operator>=(SDE other) const {return state.time >= other.state.time;}
        //bool operator!=(SDE other) const {return !(*this == other);}
        SDE_state<d> * operator*() {return &state;}
        
        using difference_type = double;
        using value_type = SDE_state<d> *;
        using pointer = const SDE_state<d> **;
        using reference = const SDE_state<d> *;
        using iterator_category = std::output_iterator_tag;
    };

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include "RNG.h"
template<const unsigned int d> struct SDE_state {
    ffloat cur[d];
    ffloat prev[d];
    ffloat time;
};

template<const unsigned int d> class SDE {
        //using SchemeEvoPolicy::evolution;
    protected:
        SDE_state<d> state;
        RNG * rng;
    public:
        SDE(const SDE_state<d>& init_cond, RNG*init_rng):state(*init_cond),rng(init_rng){}
        //SDE & operator++() { evolution(this->state,(rng->buf_cur+=d),step_width(this->state)); return *this;}
        SDE& operator++(int) {SDE& retval = *this; ++(*this); return retval;}
        //bool operator==(SDE other) const {return *this == other;}
        bool operator>=(SDE& other) const {return state.time >= other.state.time;}
        //bool operator!=(SDE other) const {return !(*this == other);}
        SDE_state<d> * operator*() {return &state;}
        using difference_type = double;
        using value_type = SDE_state<d> *;
        using pointer = const SDE_state<d> **;
        using reference = const SDE_state<d> *;
        using iterator_category = std::output_iterator_tag;
};

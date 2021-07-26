/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include<memory>
#include <list>
#include"SDE.h"
#include"Types.h"
#include"HDistribution.h"
typedef struct {
    HParams h_params;
    ffloat step_width_bound;
    ffloat barrier;
} AdaptiveHParams;
/*class HEvo {
    AdaptiveHParams params;
    protected:
        HEvo(AdaptiveHParams * params);
        void evolution(SDE_state<2> * current,ffloat * rands,double t);
};*/
class Adaptive{
    protected:
        ffloat step_width(SDE_state<2> * current); //TODO
};
class NonAdaptive{
    const ffloat step_size;
    NonAdaptive(ffloat fixed_step_size): step_size(fixed_step_size){}
    protected:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
        ffloat step_width(SDE_state<2> * current){return step_size;}
#pragma GCC diagnostic pop
};
template<typename SchemeParams,class Scheme> class HQEAnderson : SDE<2>, private Scheme{
    ffloat log_X;
    const HParams params;        
    using Scheme::step_width;
    HQEAnderson(const HParams& params_init,SchemeParams scheme_params_init, RNG*rng, const SDE_state<2> init_cond) : Scheme(scheme_params_init),SDE<2>(rng,init_cond),params(params_init),log_X(std::log(init_cond.cur[0])){}
    SDE<2>& operator++();
};

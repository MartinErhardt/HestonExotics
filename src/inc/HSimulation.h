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
    ffloat step_witdh_bound;
    ffloat barrier;
} AdaptiveHParams;
class HEvo {
    AdaptiveHParams params;
    protected:
        HEvo(AdaptiveHParams * params);
        void evolution(SDE_state<2> * current,ffloat * rands,double t);
};
class AdaptiveHEvo : HestonEvo{
    protected:
        double step_width(SDE_state<2> * current);
};
template<typename SchemeParams,class Scheme> class HMonteCarlo : SDE<2,SchemeParams, Scheme>{
    //static void calibrate(HestonParams * to_calibrate,const std::list<option>& to_calibrate_against);
    HMonteCarlo(HParams * params);
};

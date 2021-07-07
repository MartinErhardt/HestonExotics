/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include<memory>
#include <list>
#include"SDE.h"
#include"Types.h"
#include"HDistr.h"
typedef struct {
    HParams h_params;
    ffloat step_witdh_bound;
    ffloat barrier;
} AdaptiveHParams;
class HestonEvo {
    AdaptiveHParams params;
    protected:
        HestonEvo(AdaptiveHParams * params);
        void evolution(SDE_state<2> * current,ffloat * rands,double t);
};
class AdaptiveHestonEvo : HestonEvo{
    protected:
        double step_width(SDE_state<2> * current);
};
std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<option>& market_data);
template<typename SchemeParams,class Scheme> class HMC : SDE<2,SchemeParams, Scheme>{
    //static void calibrate(HestonParams * to_calibrate,const std::list<option>& to_calibrate_against);
    HMC(HParams * params);
};

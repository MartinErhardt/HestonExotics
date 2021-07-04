/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#ifndef HM_H
#define HM_H

#include"SDE.h"

typedef struct {
    double mu;
    double theta;
    double kappa;
    double xi;
    double rho;
} HestonParams;
typedef struct {
    HestonParams h_params;
    double step_witdh_bound;
    double barrier;
} AdaptiveHestonParams;
typedef struct {
    double * prices;
    double * strikes;
    //double * barriers;
} MarketData;
class HestonEvo {
    AdaptiveHestonParams params;
    protected:
        HestonEvo(AdaptiveHestonParams * params);
        void evolution(SDE_state<2> * current,double * rands,double t);
};
class AdaptiveHestonEvo : HestonEvo{
    protected:
        double step_width(SDE_state<2> * current);
};
template<typename SchemeParams,class Scheme> class HestonModel : SDE<2,SchemeParams, Scheme>{
    static void calibrate(HestonParams * to_calibrate,const std::list<option>& to_calibrate_against);
    HestonModel(HestonParams * params);
};
#endif

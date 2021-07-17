/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HCalibration.h"
#include"BSM.h"
#include"SWIFT.h"
#include<iostream>
/*typedef struct {
    ffloat mu;
    ffloat theta;
    ffloat kappa;
    ffloat xi;
    ffloat rho;
} HestonParams;
*/

std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data){
    auto to_calib=std::make_unique<HParams>();
    ffloat avg_iv=avg_imp_vol(S, market_data);
    // We start in a Black-Scholes Model see p.106 lecture notes
    *to_calib={avg_iv*avg_iv/trading_days, avg_iv*avg_iv/trading_days, 0.04,1.0,1.0};
    for(std::list<options_chain>::const_iterator opts = market_data.begin(); opts != market_data.end(); opts++){
        if(!opts->days_to_expiry) continue;
        if(opts->options.size()){
            auto current_distribution=std::make_unique<HDistribution>(*to_calib);
            auto current_SWIFT=std::make_unique<SWIFT>(*current_distribution,*SWIFT::get_parameters(*current_distribution,S,*opts));
        }
    }
    return to_calib;
}

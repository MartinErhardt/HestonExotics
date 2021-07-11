/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HCalibration.h"
#include"BSM.h"
#include"SWIFT.h"
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
        auto my_distr=std::make_unique<HDistribution>(*to_calib);
        auto current_SWIFT=std::make_unique<SWIFT>(*to_calib,S,*opts, SWIFT::get_m(*my_distr,static_cast<ffloat>(opts->days_to_expiry)/trading_days));
        current_SWIFT->get_FFT_coeffs();
    }
    return to_calib;
}

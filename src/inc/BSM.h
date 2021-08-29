/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"
#include <list>
#include<math.h>
#define EXP_LB .1
const double trading_days=static_cast<double>(5*static_cast<unsigned int>(365/7)+365%7);    //<trading days in one year
const ffloat yearly_risk_free=0.0005;                                                       //<nominal risk free rate (p.a.) 3 month T-bill yields 0.0005

/**
 * @brief computes a weighted average of the volatilty implied by the prices  in all_chains and the stock price
 * @param S stock price
 * @param all_chains list of options_chain containing European calls
 * @return weighted average of implied volatilty
 */
ffloat avg_imp_vol(const ffloat S, const std::list<options_chain>& all_chains);
/**
 * @brief computes the volatilty implied by a European call option and the stock price
 * @param S stock price
 * @param opt European call option
 * @param expi time to expiry date (in trading years)
 * @return weighted average of implied volatilty
 */
ffloat imp_vol(const ffloat S, const option& opt,ffloat expi);
ffloat call_price(const ffloat S, const ffloat K, const ffloat sigma, const double T);

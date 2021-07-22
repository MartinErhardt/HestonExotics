/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include"Types.h"
#include <list>
#include<math.h>
#define EXP_LB 0.
const double trading_days=static_cast<double>(5*static_cast<unsigned int>(365/7)+365%7);
const ffloat yearly_risk_free=0.02; //nominal risk free rate (p.a.) 3 month T-bill yields 0.0005
const ffloat risk_free= pow(1.+yearly_risk_free,1./trading_days)-1.;
ffloat avg_imp_vol(const ffloat S, const std::list<options_chain>& all_chains);

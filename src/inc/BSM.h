/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#ifndef BSM_H
#define BSM_H
#include"Types.h"
#define RISK_FREE pow(1.0005,4)-1.0 //nominal risk free rate (p.a.) 3 month T-bill yields 0.0005
ffloat imp_vol(const ffloat S, const option& opt, const ffloat precision);
#endif

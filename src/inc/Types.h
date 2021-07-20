/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include<stdint.h>
#include<vector>
#include<cstdlib>
#include <limits>

#define PRECISION 0.0000001

typedef double ffloat;
typedef struct{
    //unsigned int days_to_expiry;
    ffloat price;
    //ffloat ask;
    //ffloat bid;
    ffloat strike;
    int64_t volume;
} option;
typedef struct OptionsChain{
    std::vector<option> options;
    unsigned int days_to_expiry;
    ffloat time_to_expiry;
    ffloat max_strike;
    ffloat min_strike;
    OptionsChain(unsigned int days_until,ffloat time_until){
        options=std::vector<option>();
        min_strike=std::numeric_limits<ffloat>::max();
        max_strike=std::numeric_limits<ffloat>::lowest();
        days_to_expiry=days_until;
        time_to_expiry=time_until;
    }
}options_chain;

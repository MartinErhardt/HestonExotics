/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include<stdint.h>
#include<vector>
#include<cstdlib>
#include <limits>

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
    ffloat max_strike;
    ffloat min_strike;
    OptionsChain(unsigned int days_until){
        options=std::vector<option>();
        min_strike=std::numeric_limits<ffloat>::lowest();
        max_strike=std::numeric_limits<ffloat>::max();
        days_to_expiry=days_until;
    }
}options_chain;

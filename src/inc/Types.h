/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#include<stdint.h>
#include<vector>
#include<cstdlib>

typedef double ffloat;
typedef struct{
    unsigned int days_to_expiry;
    ffloat price;
    ffloat ask;
    ffloat bid;
    ffloat strike;
    int64_t volume;
} option;
typedef struct OptionsChain{
    std::vector<ffloat> * strikes;
    std::vector<unsigned int> * volumes;
    unsigned int days_to_expiry;
    ffloat max_strike;
    ffloat min_strike;
    ~OptionsChain(){
        free(strikes);
        free(volumes);
    }
}options_chain;

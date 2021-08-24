/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
/**
 * @file Types.h contains certain data types which are passed through different namespaces and submodules of the program.
 */
#include<stdint.h>
#include<vector>
#include<cstdlib>
#include <limits>
#include <iostream>

#if FP_SIZE==8
typedef double ffloat;  ///<financial float for more portable code depending on FP_SIZE. Currently only double is supported.
typedef uint64_t fuint; ///<unsigned integer of same size as ffloat!!
#else
#error currently unsupported
#endif
static_assert(sizeof(ffloat) == sizeof(fuint), "ffloat and fuint not of same size!");
/**
 * @brief struct encapsulating all parameters of a (call) option
 * 
 * Currently only call options are supported as there is no inherent information advantage in using put options (call-put option), which would lead to overfitting. However in cases where put options are more liquid than calls it might be a good idea.
 */
typedef struct{
    //unsigned int days_to_expiry;
    ffloat price;       //<ask price of option
    //ffloat ask;
    ffloat bid;         //<bid price of option
    ffloat strike;      //<call strike
    int64_t volume;     //<depending on program arguments either open interest of the option(contracts outstanding) or trading volume on the last(current) day. 
} option;
/**
 * @brief struct encapsulating all parameters of call option to the same expiry date and underlying, which are called a options chain.
 */
typedef struct OptionsChain{
    std::vector<option> *options;
    unsigned int days_to_expiry;    //< days to expiry
    ffloat time_to_expiry;          //<time_to_expiry (in years)
    ffloat max_strike;              //< minimum strike of all options
    ffloat min_strike;              //< maximum strike of all options
    OptionsChain(unsigned int days_until,ffloat time_until){
        options=new std::vector<option>();
        min_strike=std::numeric_limits<ffloat>::max();
        max_strike=std::numeric_limits<ffloat>::lowest();
        days_to_expiry=days_until;
        time_to_expiry=time_until;
    }
    ~OptionsChain(){delete options;}
}options_chain;

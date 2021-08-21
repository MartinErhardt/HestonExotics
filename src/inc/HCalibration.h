/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include"HDistribution.h"
#include<memory>
#include<list>
#include"SWIFT.h"
#include"BSM.h"
typedef struct ED:traced<ED>{
    const options_chain& opts;
    HDistribution* distr;
    SWIFT* pricing_method;
    //unsigned int underpriced;
    ED(const options_chain& opts,HDistribution* init_distr, SWIFT* init_pricing_method): opts(opts),distr(init_distr),pricing_method(init_pricing_method){}
    //ED(const options_chain& opts,HParams p, ffloat expi): opts(opts),distr(new HDistribution(p,expi)){}
    ~ED(){delete distr; delete pricing_method;}
    //ED(const ED& to_copy){
    //    std::copy
    //}
} expiry_data;
typedef struct AS{
    ffloat S;
    ffloat* real_prices;
    std::list<expiry_data>& exp_list;
    ~AS(){delete &exp_list;} //moved when push_back 
    
} adata_s;
std::ostream& operator<<(std::ostream& out, adata_s const& as);
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata);
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata);

std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data);

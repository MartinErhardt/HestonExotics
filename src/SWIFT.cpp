/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"SWIFT.h"
#include"BSM.h"
#include <iostream>


SWIFT::SWIFT(const HDistribution& distr, const ffloat stock_price, const options_chain& opts,const unsigned int m_d) : distr(distr), S(stock_price), m(m_d){
    ffloat max=risk_free*opts.days_to_expiry+std::log(stock_price/opts.min_strike);
    ffloat min=risk_free*opts.days_to_expiry+std::log(stock_price/opts.max_strike);
    
    ffloat c=std::abs(distr.first_order_moment(opts.days_to_expiry))+10.*std::sqrt(std::fabs(distr.second_order_moment(opts.days_to_expiry))+std::sqrt(std::abs(distr.fourth_order_moment(opts.days_to_expiry))));
    ffloat exp2_m=std::exp2(m);
    this->k_1=ceil(exp2_m*(min - c));
    this->k_2=floor(exp2_m*(max + c));
    ffloat iota_density=ceil(std::log2(M_PI*std::abs(k_1-k_2)))-1;
    this->J=std::exp2(iota_density-1);
    
    //int cas;
    //for(i=1;(i<=k_2-k_1)||(cas=(i<=-k_1));i++) payoffs[-k_1+cas*i-(1-cas)*(i+k_1))]=(out[i]+(cas*2-1)*out2[i-1])*0.5;
    
}
unsigned int SWIFT::get_m(const HDistribution& distr,ffloat tau){
    //p=params;
    // find m by bisection
    unsigned int i_val=1;
    unsigned int l_val;
    unsigned int u_val;
    unsigned int m_val;
    unsigned int cur_error;
    //std::cout<<"start_val: "<<i_val<<'\n';
    if (distr.int_error(i_val,tau)<PRECISION){
        l_val=i_val;
        //std::cout<<"start_val: "<<i_val<<'\n';
        do{i_val*=2;
            std::cout<<"i_val: "<<i_val<<'\n';
            cur_error=distr.int_error(i_val,tau);
            std::cout<<"i_val_after: "<<i_val<<'\n';
        }
        while (cur_error<PRECISION && i_val);
        //std::cout<<"end_val: "<<i_val<<'\n';
        u_val=i_val;
    } else{
        u_val=i_val;
        do{ i_val=static_cast<unsigned int>(i_val*.5);
            //std::cout<<"i_val2: "<<i_val<<'\n';
            
        }
        while(distr.int_error(i_val,tau)>=PRECISION && (i_val>0));
        l_val=i_val;
    }
    while(u_val-l_val>1){
        m_val=static_cast<unsigned int>((l_val+u_val)/2);
        cur_error=distr.int_error(i_val,tau);
        if(cur_error<PRECISION) l_val=m_val;
        else u_val=m_val;
    }
    return u_val;
}
void SWIFT::get_payoff(){
    
}

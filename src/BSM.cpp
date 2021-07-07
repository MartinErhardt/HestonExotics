/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"BSM.h"
#include<math.h>
#include <list>
#include <iostream>

#define PRECISION 0.0000001

ffloat norm_cdf(const ffloat value);
ffloat d_j(const int j, const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T);
ffloat call_price(const ffloat S, const ffloat K, const ffloat sigma, const double T);
ffloat imp_vol(const ffloat S, const option& opt,unsigned int days_to_expiry);

ffloat norm_cdf(ffloat value){
   return 0.5 * erfc(-value * M_SQRT1_2);
}
ffloat d_j(const int j, const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T){
  return (log(S/K)+(r+(pow(-1,j-1))*0.5*sigma*sigma)*T)/(sigma*(pow(T,0.5)));
}
ffloat call_price(const ffloat S, const ffloat K, const ffloat sigma, const double T){
    return S*norm_cdf(d_j(1,S,K,risk_free,sigma,T))-K*exp(-risk_free*T) * norm_cdf(d_j(2,S,K,risk_free,sigma, T));
}
ffloat imp_vol(const ffloat S, const option& opt,unsigned int days_to_expiry){
    ffloat i_val=(1/64.0);
    ffloat l_val;
    ffloat u_val;
    ffloat m_val;
    ffloat p;
    const double expi=static_cast<double>(days_to_expiry);
    double lower_bound=S-std::exp(-risk_free*expi)*opt.strike;
    if(opt.price<=lower_bound) return -1.0;
    if(opt.price>=S) return -2.0;
    if (call_price(S,opt.strike,i_val, expi)<opt.price){
        l_val=i_val;
        do i_val*=2;
        while (call_price(S,opt.strike, i_val, expi)<opt.price&& i_val>0);
        while(i_val==0);
        u_val=i_val;
    } else{
        u_val=i_val;
        do i_val*=0.5;
        while (call_price(S,opt.strike, i_val, expi)>=opt.price && i_val >0);
        while(i_val==0);
        l_val=i_val;
    }
    while(std::abs(u_val-l_val)>=PRECISION){
        m_val=(l_val+u_val)/2;
        p=call_price(S,opt.strike,m_val, expi);
        if(p<opt.price) l_val=m_val;
        else u_val=m_val;
    }
    return m_val;
}
ffloat avg_imp_vol(const ffloat S, const std::list<options_chain>& all_chains){
    ffloat avg_vol=0;
    ffloat denom=0;
    double imp_v;
    for(std::list<options_chain>::const_iterator cur_chain = all_chains.begin(); cur_chain != all_chains.end(); cur_chain++){
        for(const option& opt: cur_chain->options){
            if((imp_v=imp_vol(S,opt,cur_chain->days_to_expiry))>=0){
                avg_vol+=opt.volume*imp_v;
                denom+=opt.volume;
            } else std::cout<<"Can't calculate implied volatility (price below lower arbitrage-bound)\n";
            std::cout<<"stock price: "<<S<<"\tprice: "<<opt.price<<"\tstrike: "<<opt.strike<<"\trisk free rate: "<<risk_free<<"\ttrading volume: "<<opt.volume<<"\timplied volatility: "<<imp_v<<'\n';
        }
        
    }
    return avg_vol/denom;
}

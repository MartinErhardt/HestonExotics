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
ffloat imp_vol(const ffloat S, const option& opt);

ffloat norm_cdf(ffloat value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}
ffloat d_j(const int j, const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T) 
{
  return (log(S/K)+(r+(pow(-1,j-1))*0.5*sigma*sigma)*T)/(sigma*(pow(T,0.5)));
}
ffloat call_price(const ffloat S, const ffloat K, const ffloat sigma, const double T) 
{
    return S*norm_cdf(d_j(1,S,K,risk_free,sigma,T))-K*exp(-risk_free*T) * norm_cdf(d_j(2,S,K,risk_free,sigma, T));
}
ffloat imp_vol(const ffloat S, const option& opt)
{
    ffloat i_val=(1/64.0);
    ffloat l_val;
    ffloat u_val;
    ffloat m_val;
    ffloat p;
    const double expi=static_cast<double>(opt.days_to_expiry);
    double lower_bound=S-std::exp(-risk_free*expi)*opt.strike;
    if(opt.price<=lower_bound&&lower_bound>0) return -1.0; //TODO Check upper bound
    if(opt.price>=S) return -2.0;
    if (call_price(S,opt.strike,i_val, expi)<opt.price){
        l_val=i_val;
        do i_val*=2;
        while (call_price(S,opt.strike, i_val, expi)<opt.price);
        u_val=i_val;
    } else{
        u_val=i_val;
        do i_val*=0.5;
        while (call_price(S,opt.strike, i_val, expi)>=opt.price);
        l_val=i_val;
    }
    while(std::abs(u_val-l_val)>=PRECISION){
        //std::cout<<"u_val: "<<u_val<<"\t l_val: "<<l_val<<"\tdiff:"<<std::abs(u_val-l_val)<<'\n';
        m_val=(l_val+u_val)/2;
        p=call_price(S,opt.strike,m_val, expi);
        if(p<opt.price) l_val=m_val;
        else u_val=m_val;
    }
    return m_val;
}
ffloat avg_imp_vol(const ffloat S, const std::list<option>& opts){
    ffloat avg_vol=0;
    ffloat denom=0;
    double imp_v;
    for(std::list<option>::const_iterator opt = opts.begin(); opt != opts.end(); opt++){
        if((imp_v=imp_vol(S,*opt))>=0){
            avg_vol+=opt->volume*imp_vol(S,*opt);
            denom+=opt->volume;
        } else std::cout<<"Can't calculate implied volatility (price below lower arbitrage-bound)\n";
        std::cout<<"stock price: "<<S<<"\tprice: "<<opt->price<<"\tstrike: "<<opt->strike<<"\trisk free rate: "<<risk_free<<"\ttrading volume: "<<opt->volume<<"\tdays to expiry: "<<opt->days_to_expiry<<"\timplied volatility: "<<imp_v<<"\t spread: "<<opt->ask-opt->bid<<'\n';
    }
    return avg_vol/denom;
}

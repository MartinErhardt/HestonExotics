/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"BSM.h"
#include<math.h>

ffloat norm_cdf(const ffloat value);
ffloat d_j(const int j, const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T);
ffloat call_price(const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T);
ffloat norm_cdf(ffloat value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}
ffloat d_j(const int j, const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T) 
{
  return (log(S/K) + (r + (pow(-1,j-1))*0.5*sigma*sigma)*T)/(sigma*(pow(T,0.5)));
}
ffloat call_price(const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T) 
{
  return S * norm_cdf(d_j(1, S, K, r, sigma, T))-K*exp(-r*T) * norm_cdf(d_j(2, S, K, r, sigma, T));
}
ffloat imp_vol(const ffloat S, const option& opt, const ffloat precision)
{
    ffloat i_val=0.125;
    ffloat l_val;
    ffloat u_val;
    ffloat m_val;
    ffloat p;
    const double expi=opt.days_to_expiry/365;
    if (call_price(S,opt.strike,RISK_FREE,i_val, opt.days_to_expiry)<opt.price){
        l_val=i_val;
        while (call_price(S,opt.strike,RISK_FREE,(i_val*=2), expi)<opt.price);
        u_val=i_val;
    } else{
        u_val=i_val;
        while (call_price(S,opt.strike,RISK_FREE,(i_val*=0.5), expi)>=opt.price);
        l_val=i_val;
    }
    while(u_val-l_val>=precision){
        m_val=(l_val+u_val)/2;
        p=call_price(S,opt.strike,RISK_FREE,m_val, expi);
        if(p<opt.price) l_val=m_val;
        else u_val=m_val;
    }
    return m_val;
}

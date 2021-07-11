#pragma once
#include"HDistribution.h"

class SWIFT{
    const HDistribution& distr;
    const ffloat S;
    const unsigned int m;
    ffloat k_1;
    ffloat k_2;
    ffloat J;
public:
    SWIFT(const HDistribution& new_distr, const ffloat stock_price, const options_chain& opts,const unsigned int m_d);
    std::vector<ffloat> price_opts(ffloat S, std::vector<double> const& Ks);
    std::vector<std::vector<ffloat>> price_opts_grad(std::vector<double> const& Ks);
    void get_payoff();
    static unsigned int get_m(const HDistribution& distr,ffloat tau);
};

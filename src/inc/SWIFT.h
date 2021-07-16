#pragma once
#include"HDistribution.h"

class SWIFT{
    const HDistribution& distr;
    const ffloat S;
    const unsigned int m;
    int exp2_m;
    ffloat sqrt_exp2_m;
    ffloat lower;
    ffloat upper;
    int k_1;
    int k_2;
    int J;
    std::vector<std::complex<ffloat>>*density_coeffs;
    std::complex<ffloat> H(ffloat y, int j);
public:
    SWIFT(const HDistribution& new_distr, const ffloat stock_price, const options_chain& opts,const unsigned int m_d);
    std::vector<ffloat> price_opts(ffloat S, std::vector<double> const& Ks);
    std::vector<std::vector<ffloat>> price_opts_grad(std::vector<double> const& Ks);
    void get_FFT_coeffs();
    static unsigned int get_m(const HDistribution& distr,ffloat tau);
    ~SWIFT();
};

#pragma once
#include"HDistribution.h"
#include<memory>
#include<list>

typedef struct{
    const unsigned int m;
    const unsigned int exp2_m;
    const ffloat sqrt_exp2_m;
    const ffloat lower;
    const ffloat upper;
    const int k_1;
    const int k_2;
    const unsigned int J;
} swift_parameters;
class SWIFT{
    HDistribution& distr;
    const swift_parameters my_params;
    std::vector<std::complex<ffloat>>& density_coeffs;
    typedef struct CacheEntry{
        const options_chain& to_price;
    } cache_entry;
    std::list<cache_entry> results_cache;
    void get_FFT_coeffs();
public:
    SWIFT(HDistribution& init_distr, const swift_parameters& init_params);
    std::vector<ffloat> price_opts(ffloat S, std::vector<double> const& Ks);
    std::vector<std::vector<ffloat>> price_opts_grad(std::vector<double> const& Ks);
    static std::unique_ptr<swift_parameters> get_parameters(const HDistribution& distr,const ffloat stock_price, const options_chain& opts);
    ~SWIFT();
};

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include"HDistribution.h"
#include<memory>
#include<list>
#include <eigen3/Eigen/Dense>
typedef struct SwiftParameters{
    const unsigned int m;
    const unsigned int exp2_m;
    const ffloat sqrt_exp2_m;
    const ffloat lower;
    const ffloat upper;
    const int k_1;
    const int k_2;
    const unsigned int J;
    ffloat u(const unsigned int i) const;
} swift_parameters;
std::ostream& operator<<(std::ostream& out, swift_parameters const& sp);
#if FP_SIZE==8
typedef Eigen::MatrixXcd MatC;
#else
#error no Matrix format found
#endif
class SWIFT:traced<SWIFT>{
public:
    const swift_parameters my_params;
private:
    std::vector<std::complex<ffloat>>& density_coeffs;
    typedef struct CacheEntry{
        MatC results;
        const options_chain& to_price;
        CacheEntry(const HDistribution& distr,const SWIFT& swift_obj, const options_chain& to_price_init,const ffloat stock_price);
        //~CacheEntry();
    } cache_entry;
    std::list<cache_entry>& results_cache;
    SWIFT::cache_entry * get_precached(const HDistribution& distr,const ffloat S, const options_chain& opts);
    //Eigen::MatrixXd * pricing_matrix = nullptr;
    void get_FFT_coeffs();
public:
    SWIFT(const swift_parameters& init_params);
    void flush_cache();
    void price_opts(const HDistribution& distr,const ffloat S, const options_chain& opts,ffloat** out_array,ffloat*end);
    void price_opts_grad(const HDistribution& distr,const ffloat S, const options_chain& opts, ffloat** out_array,ffloat*end);
    static std::unique_ptr<swift_parameters> get_parameters(const HDistribution& distr,const ffloat stock_price, const options_chain& opts);
    ~SWIFT();
};

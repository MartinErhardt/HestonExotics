/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include"HDistribution.h"
#include<memory>
#include<list>
#include <eigen3/Eigen/Dense>
/** @file SWIFT.h contains an implementation of the SWIFT pricing method, which differs from the one provided by Eudald Romo Grau in (1) the following ways
 * 1. Iterative method as described in (2) being applied to determine the wavelet scale parameter m
 * 2. FFTW is applied to compute both Fourier transformations including the second one described in the paper.
 * 3. the last sum computation in (3) for price and all gradient components is simultaneously performed as Matrix-Vector multiplication using Eigen, with the result being cached
 * 4. Simplification of Class structure
 * @see https://github.com/eudaldrg
 * @see https://papers.ssrn.com/sol3/Delivery.cfm/SSRN_ID2823900_code2491804.pdf?abstractid=2705699&mirid=1
 * @see https://arxiv.org/pdf/2103.01570.pdf
 */
/**
 * @brief Contains various relevant parameters to the SWIFT algorithm
 * 
 * Very closely resembles SwiftParameters in Eudald Romo Graus implentation with the difference being the iterative method as described in (1) being applied to determine the wavelet scale parameter m.
 * 
 * @see https://papers.ssrn.com/sol3/Delivery.cfm/SSRN_ID2823900_code2491804.pdf?abstractid=2705699&mirid=1
 * @see https://papers.ssrn.com/sol3/Delivery.cfm/SSRN_ID2586531_code1988671.pdf?abstractid=2585529&mirid=1
 * @see https://github.com/eudaldrg/SWIFTOptionCalibration/blob/master/SWIFT/swift_parameters.h
 */
typedef struct SwiftParameters{
    unsigned int m;         //<wavelet scale parameter m
    unsigned int exp2_m;    //<two to the power of m
    ffloat sqrt_exp2_m;     //<root of two to the power of m
    ffloat lower;           //<lower integral bound
    ffloat upper;           //<upper integral bound
    int k_1;                //<lower sum index
    int k_2;                //<upper sum index
    unsigned int J;         //<Fourier size(both density and payoff coefficients)
    /**
     * this struct method computes the frequently used index u_i
     * @param i positive sum index
     * @return pi/(2*J)*(2*i+1) (Note that we index starting in 0, while the paper indexes starting in 1)
     */
    ffloat u(const unsigned int i) const;
    /**
     * Constructor dynamically computing the coefficients
     *
     * @param distr Heston model distribution
     * @param stock_price stock price
     * @param opts options to price
     */
    SwiftParameters(const HDistribution& distr,const ffloat stock_price, const options_chain& opts);
    /**
     * Simple initializer constructor
     * 
     * This constructor is used in UnitTest.cpp
     * 
     * @param m wavelet scale parameter m
     * @param exp2_m two to the power of m
     * @param sqrt_exp2_m root of two to the power of m
     * @param lower lower integral bound
     * @param upper upper integral bound
     * @param k_1 lower sum index
     * @param k_2 upper sum index
     * @param J Fourier size(both density and payoff coefficients)
     */
    SwiftParameters(const unsigned int m,unsigned int exp2_m,ffloat sqrt_exp2_m,ffloat lower,ffloat upper,int k_1,int k_2,unsigned int J):m(m),exp2_m(exp2_m),sqrt_exp2_m(sqrt_exp2_m),lower(lower),upper(upper),k_1(k_1),k_2(k_2),J(J){}
} swift_parameters;
/**
 * custom shift operator used for debugging in the standard output, output format:
 * m: (m)[tab]exp2_m:(exp2_m)[tab]sqrt_exp2_m: (sqrt_exp2_m)[tab]lower: (lower)[tab]upper: (upper)[tab]k_1: (k_1)[tab]k_2: (k_2)[tab]J: (J)
 */
std::ostream& operator<<(std::ostream& out, swift_parameters const& sp);
#if FP_SIZE==8
typedef Eigen::MatrixXcd MatC; //<typedef for Matrixformat here MatrixXcd, because only ffloat double is supported
#else
#error no Matrix format found
#endif
/**
 * @brief Given constant swift parameters this class encapsulates the density coefficients, which are calculated only once using FFT
 * 
 * For each options chain(all options to same expiry date) the results of both price and gradient information are stored in a results cache, such that the time intensive calculation of the characteristic fanction has to be done only once during each Levenberg-Marquardt iteration.
 * 
 * This class most closely resembles the class SwiftInvariantData in the implementation provided by Eudald Romo Grau, but also includes functionality of FFTPayoffCalculator, which do not seem to offer any advantage as an abstraction. It also includes functionality of the classes OptionContract and EuropeanOptionContract, though this might change later.
 * 
 * @see https://github.com/eudaldrg/SWIFTOptionCalibration/blob/master/SWIFT/quick_callibration_swift.h 
 * @see https://github.com/eudaldrg/SWIFTOptionCalibration/blob/master/SWIFT/payoff_coefficients_calculators.cpp
 * @see https://github.com/eudaldrg/SWIFTOptionCalibration/blob/master/SWIFT/option_contracts.h
 */
class SWIFT{
public:
    const swift_parameters my_params; //< constant parameters
private:
    std::vector<std::complex<ffloat>> density_coeffs; //<note that contrary to the usual intuition with references this reference is owned by the instance of this class itself not some other part of the program. This is done to allow usage of the overloaded [] operator.
    /**
     * @brief Private CacheEntry datatype
     * Stores Pricing and Gradient information in 6xto_price.options Eigen matrix
     */
    typedef struct CacheEntry{
        MatC results; //<6xto_price.options Eigen matrix
        const options_chain& to_price; //< Options chain to be priced
        /**
         * Constructor for cache entry
         * 
         * @param distr Heston model distribution
         * @param swift_obj instance of SWIFT used to do the calculation. This always refers to the same instance in which the results_cache in which this cache_entry is stored is part of. 
         * @param to_price_init options chain the pricing and gradient information of which is to be determined.
         * @param stock_price stock price
         */
        CacheEntry(const HDistribution& distr,const SWIFT& swift_obj, const options_chain& to_price_init,const ffloat stock_price);
        //~CacheEntry();
    } cache_entry;
    std::list<cache_entry> results_cache;//<note that contrary to the usual intuition with references this reference is owned by the instance of this class itself not some other part of the program.
    /**
     * private method returning a pointer to the cache_entry in which pricing and gradient information is stored if available and returns a null pointer otherwise
     * 
     * @param distr Heston model distribution to which we want to determine pricing or gradient information
     * @param S price of stock
     * @param opts options chain to price/find gradient
     * @return C-style pointer to cache_entry
     */
    SWIFT::cache_entry * get_precached(const HDistribution& distr,const ffloat S, const options_chain& opts);
public:
    /**
     * Constructor computing the density coefficients & Initizialization of const swift_parameters through ctor
     * 
     * Unlike the original implementation by Eudald Romo Grau this constructor computes both Fourier transformations including the last one efficiently using FFTW
     * 
     * @param init_params to be copied into my_params
     */
    SWIFT(const swift_parameters& init_params);
    /**
     * flush the internal cache
     * 
     * do this, when the distribution(more specifically its parameters) to which price_opts will be applied have changed. Not necessary though.
     * 
     * @sideeffect empty clear results_cache
     */
    void flush_cache();
    /**
     * This method looks for a suitable cache entry to obtain pricing information. Otherwise it generates one.
     * 
     * Note that the C-style call syntax of this method stems from its' application to C-library dlevmar
     * 
     * @param distr Distribution containing the Heston Parameters and expiry information
     * @param S price of the underlying
     * @param opts option chain to price
     * @param out_array start of the output array
     * @param end end of allocated output buffer 
     * @throws std::runtime_exception if the buffer in *out_array is too small to accomadate all of the pricing information
     * @sideeffect **out_array is filled with the prices of the options in opt.options as determined by the SWIFT method in exactly the order to be found there. Furthermore *out_array is updated to point at the first unwritten entry before end.
     */
    void price_opts(const HDistribution& distr,const ffloat S, const options_chain& opts,ffloat** out_array,ffloat*end);
    /**
     * This method looks for a suitable cache entry to obtain gradient information. Otherwise it generates one.
     * 
     * Note that the C-style call syntax of this method stems from its' application to the C-library dlevmar
     * 
     * @param distr Distribution containing the Heston Parameters and expiry information
     * @param S price of the underlying
     * @param opts option chain to obtain gradient information on.
     * @param out_array start of the output array
     * @param end end of allocated output buffer 
     * @throws std::runtime_exception if the buffer in *out_array is too small to accomodate all of the gradient information
     * @sideeffect out_array is filled with the gradient of the options in opt.options as determined by the SWIFT method in exactly the order to be found there. Furthermore *out_array is updated to point at the first unwritten entry before end. The various components are ordered as in HParams.
     */
    void price_opts_grad(const HDistribution& distr,const ffloat S, const options_chain& opts, ffloat** out_array,ffloat*end);
};

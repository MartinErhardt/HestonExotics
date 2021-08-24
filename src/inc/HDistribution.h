#pragma once 
#include"Types.h"
#include<complex>
#include<vector>
#include<iostream>
/**
 * @brief struct containing all parameters to calibrate for
 */
typedef struct HP{
    ffloat v_0;     //< inital varince
    ffloat v_m;     //< long term variance
    ffloat rho;     //< correlation between spot and volatility
    ffloat kappa;   //< mean reversion rate
    ffloat sigma;   //< variance of volatility
    /**
     * Compares the content of parameters l and r
     * @param l left-side Heston parameter
     * @param r right-side Heston parameter
     * @return equality in all components.
     */
    friend bool operator==(const HP& l, const HP& r){
        return (l.v_0==r.v_0) && (l.v_m==r.v_m) && (l.rho==r.rho) && (l.kappa==r.kappa) && (l.sigma==r.sigma);
    };
} HParams;
/**
 * @brief Models the probability distribution of the price process of the underlying depending on HParams p, time tau outstanding and the risk_free rate r
 * 
 * The methods of this class compute the characteristic function, the gradient of the characteristic function the first, second and fourth order cumulants and the integration error w.r.t to the SWIFT method.
 */
class HDistribution{
    
    struct Helpers{
        /**
        * struct containing various Helper variables, the first seven of which are described in https://arxiv.org/pdf/2103.01570.pdf . The others are precomputed as they take up a certain amount of time to compute.
        */
        std::complex<ffloat> xi;
        std::complex<ffloat> d;
        std::complex<ffloat> A_1;
        std::complex<ffloat> A_2;
        std::complex<ffloat> A;
        std::complex<ffloat> B;
        std::complex<ffloat> D;
        
        std::complex<ffloat> sinh_v;        //<sinus hyperbolicus
        std::complex<ffloat> cosh_v;        //<cosinus hyperbolicus
        std::complex<ffloat> exp_kappa_tau; //<$e^{\frac{kappa tau}{2}}$
        Helpers(const HParams& p,const std::complex<ffloat> u,const ffloat tau);
    };
    typedef Helpers helpers;
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u,const helpers& hlp,std::complex<ffloat>chf_val) const;
    std::vector<std::complex<ffloat>> chf_chf_grad(const std::complex<ffloat> ,const helpers& hlp,std::complex<ffloat>chf_val) const;
    std::complex<ffloat> chf(const std::complex<ffloat> u,const helpers& hlp) const;
public:
    HParams p;              //< Heston model parameters
    const ffloat tau;       //< time at which the stochastic process is observed(usually time to expiry in years)
    const ffloat risk_free; //< risk free rate
    HDistribution(HParams params,const ffloat init_tau,const ffloat risk_free_init): p(HParams(params)),tau(init_tau),risk_free(risk_free_init){};
    /**
     * Computes both the characteristic function and all of its partial derivatives
     * @param u
     * @return evaluation of the characteristic function in u
     */
    std::complex<ffloat> chf(const std::complex<ffloat> u) const;
    /**
     * Computes both the characteristic function and all of its partial derivatives
     * @param u
     * @return complex vector of length five w/ entries containing the evaluations of the partial derivatives in u in the order of HParams
     */
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u) const;
    /**
     * Computes both the characteristic function and all of its partial derivatives
     * @param u
     * @return complex vector of length six: The first entry contains the evaluations of the characteristic function in u, the others those of the partial derivatives in u in the order of HParams.
     */
    std::vector<std::complex<ffloat>> chf_chf_grad(const std::complex<ffloat> u) const;
    /**
     * integration error with wave length parameter m
     * 
     * @param trunc_m wavelength parameter m
     * @return error upper bound
     */
    ffloat int_error(const unsigned int trunc_m) const;
    /**
     * Computes expectation
     * 
     * @author Eudald Romo Grau (eudald.romo@xanadutrading.com)
     * @return expectation
     */
    ffloat first_order_moment() const;
    /**
     * Computes second order cumulant
     * 
     * @author Eudald Romo Grau (eudald.romo@xanadutrading.com)
     * @return distributional moment of second order
     */
    ffloat second_order_moment() const;
    /**
     * Computes fourth order cumulant
     * 
     * @author Eudald Romo Grau (eudald.romo@xanadutrading.com)
     * @return distributional moment of fourth order
     */
    ffloat fourth_order_moment() const;
    //~HDistribution(){std::cout<<"I'm being destructed:((\n";};
};


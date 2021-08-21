#pragma once 
#include"Types.h"
#include<complex>
#include<vector>
#include<iostream>

typedef struct HP{
    ffloat v_0; // inital varince
    ffloat v_m; // long term variance
    ffloat rho; // correlation between spot and volatility
    ffloat kappa; // mean reversion rate
    ffloat sigma; // variance of volatility
    friend bool operator==(const HP& l, const HP& r){
        //std::cout<<"== overload successful\n";
        return (l.v_0==r.v_0) && (l.v_m==r.v_m) && (l.rho==r.rho) && (l.kappa==r.kappa) && (l.sigma==r.sigma);
    };
} HParams;

class HDistribution:traced<HDistribution>{
    struct Helpers{
        std::complex<ffloat> xi;
        std::complex<ffloat> d;
        std::complex<ffloat> A_1;
        std::complex<ffloat> A_2;
        std::complex<ffloat> A;
        std::complex<ffloat> B;
        std::complex<ffloat> D;
        
        std::complex<ffloat> sinh_v;
        std::complex<ffloat> cosh_v;
        std::complex<ffloat> exp_kappa_tau;
        Helpers(const HParams& p,const std::complex<ffloat> u,const ffloat tau);
    };
    typedef Helpers helpers;
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u,const helpers& hlp,std::complex<ffloat>chf_val) const;
    std::vector<std::complex<ffloat>> chf_chf_grad(const std::complex<ffloat> ,const helpers& hlp,std::complex<ffloat>chf_val) const;
    std::complex<ffloat> chf(const std::complex<ffloat> u,const helpers& hlp) const;
public:
    HParams p;
    const ffloat tau;
    const ffloat risk_free;
    HDistribution(HParams params,const ffloat init_tau,const ffloat risk_free_init): p(HParams(params)),tau(init_tau),risk_free(risk_free_init){};
    std::complex<ffloat> chf(const std::complex<ffloat> u) const;
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u) const;
    std::vector<std::complex<ffloat>> chf_chf_grad(const std::complex<ffloat> u) const;
    ffloat int_error(const unsigned int trunc_m) const;
    ffloat first_order_moment() const;
    ffloat second_order_moment() const;
    ffloat fourth_order_moment() const;
    //~HDistribution(){std::cout<<"I'm being destructed:((\n";};
};


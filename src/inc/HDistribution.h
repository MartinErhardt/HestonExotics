#pragma once 
#include"Types.h"
#include<complex>
#include<vector>

typedef struct {
    ffloat v_0; // inital varince
    ffloat v_m; // long term variance
    ffloat rho; // correlation between spot and volatility
    ffloat kappa; // mean reversion rate
    ffloat sigma; // variance of volatility
} HParams;

class HDistribution{
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
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u,const ffloat tau,const helpers& hlp,std::complex<ffloat>chf_val) const;
    std::complex<ffloat> chf(const std::complex<ffloat> u,const ffloat tau,const helpers& hlp) const;
public:
    HParams p;
    HDistribution(HParams params){p=params;}
    std::complex<ffloat> chf(const std::complex<ffloat> u,const ffloat tau) const;
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u,const ffloat tau) const;
    ffloat int_error(const unsigned int trunc_m, const ffloat tau) const;
    ffloat first_order_moment(ffloat T) const;
    ffloat second_order_moment(ffloat T) const;
    ffloat fourth_order_moment(ffloat T) const;
};
void distr_test();

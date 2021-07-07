#pragma once 
#include"Types.h"
#include<complex>
#include<vector>

typedef struct {
    ffloat v_0;
    ffloat v_m;
    ffloat rho;
    ffloat kappa;
    ffloat sigma;
/*
    	double kappa;         // mean reversion rate
	double v_m;     // long term variance
	double m_sigma;     // variance of volatility
	double m_rho;       // correlation between spot and volatility
	double m_v0;        // initial variance
*/
} HParams;

class HDistribution{
    HParams p;
    ffloat T;

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
        public: Helpers(const HParams& p,const std::complex<ffloat> u,const ffloat tau);
    };
    typedef Helpers helpers;
public:
    std::complex<ffloat> chf(const std::complex<ffloat> u,const ffloat tau,const helpers& hlp);
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u,const ffloat tau,const helpers& hlp,std::complex<ffloat>chf_val);
    std::complex<ffloat> chf(const std::complex<ffloat> u,const ffloat tau);
    std::vector<std::complex<ffloat>> chf_grad(const std::complex<ffloat> u,const ffloat tau);

    /*
    std::complex<ffloat> chf2(const std::complex<ffloat> u,const ffloat tau);
    std::vector<std::complex<ffloat>> chf_grad2(const std::complex<ffloat> u,const ffloat tau);
    */
    HDistribution(HParams params);
    ffloat first_order_moment();
    ffloat second_order_moment();
    ffloat fourth_order_moment();
};

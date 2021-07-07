#pragma once 
#include"Types.h"
#include<complex>
#include<vector>

typedef struct {
    ffloat nu_i;
    ffloat nu_m;
    ffloat rho;
    ffloat kappa;
    ffloat sigma;
    
} HParams;
class HDistribution{
    HParams p;
    double expiry;
public:
    HDistribution(HParams params);
    std::complex<ffloat> chf();
    void chf_grad(std::vector<std::complex<ffloat>>& grad);
    void chf_grad(std::vector<std::complex<ffloat>>& grad,const std::complex<ffloat> chf);
    ffloat first_order_moment();
    ffloat second_order_moment();
    ffloat fourth_order_moment();
};

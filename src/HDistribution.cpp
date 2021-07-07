#include"HDistribution.h"
#include"BSM.h"
#include <cmath>
#include<complex.h>
using namespace std::complex_literals;

HDistribution::Helpers::Helpers(const HParams& p,const std::complex<ffloat> u,const ffloat tau){
        xi=p.kappa+p.sigma*p.rho*1i*u;
        std::complex<ffloat> fac=u*u+1i*u;
        d=std::sqrt(xi*xi+p.sigma*p.sigma*fac); //Why does this work?
        const std::complex<ffloat> exp_1=std::exp(d*tau*.5);
        const std::complex<ffloat> exp_2=1./exp_1;
        sinh_v=(exp_1-exp_2)*.5; //numerically instable
        cosh_v=(exp_1+exp_2)*.5;
        A_1=fac*sinh_v;
        A_2=d/p.v_0*cosh_v+xi/p.v_0*sinh_v;
        A=A_1/A_2;
        exp_kappa_tau=std::exp(p.kappa*tau*.5);
        B=d*exp_kappa_tau/(p.v_0*A_2);
        D=std::log((2.*d)/(d+xi+(d-xi)*std::exp(-d*tau)))+(p.kappa-d)*tau*.5;
}
std::complex<ffloat> HDistribution::chf(const std::complex<ffloat> u,const ffloat tau){
    HDistribution::helpers hlp(this->p,u,tau);
    return chf(u,tau,hlp);
}
std::complex<ffloat> HDistribution::chf(const std::complex<ffloat> u,const ffloat tau,const helpers& hlp){
    return std::exp(-1i*u*risk_free*tau+1i*p.kappa*p.v_m*p.rho*u/p.sigma-hlp.A+2.*hlp.D*p.kappa*p.v_m/(p.sigma*p.sigma));
}
std::vector<std::complex<ffloat>> HDistribution::chf_grad(const std::complex<ffloat> u,const ffloat tau){
    HDistribution::helpers hlp(this->p,u,tau);
    std::complex<ffloat> chf_val_arg=chf(u,tau,hlp);
    return chf_grad(u,tau,hlp,chf_val_arg);
}
std::vector<std::complex<ffloat>> HDistribution::chf_grad(const std::complex<ffloat> u,const ffloat tau,const helpers& hlp,const std::complex<ffloat>chf_val){
    std::complex<ffloat> d_rho=hlp.xi*p.sigma*1i*u/hlp.d;
    std::complex<ffloat> A_2_rho=p.sigma*1i*u*(2.+hlp.xi*tau)/(2.*hlp.d*p.v_0)*(hlp.xi*hlp.cosh_v+hlp.d*hlp.sinh_v);
    std::complex<ffloat> B_rho=hlp.exp_kappa_tau/p.v_0*(d_rho/hlp.A_2-hlp.d/(hlp.A_2*hlp.A_2)*A_2_rho);
    std::complex<ffloat> A_1_rho=(1i*u*(u*u-1i*u)*tau*hlp.xi*p.sigma)/(2.*hlp.d)*hlp.cosh_v;
    std::complex<ffloat> A_rho=A_1_rho/hlp.A_2-hlp.A/hlp.A_2*A_2_rho;
    //std::complex<ffloat> A_kappa=-1i/(p.sigma*u)*A_rho;
    std::complex<ffloat> B_kappa=-1i/(p.sigma*u)*B_rho+hlp.B*tau*.5;
    std::complex<ffloat> d_sigma=(hlp.d/p.sigma-1./hlp.xi)*d_rho+p.sigma*u*u/hlp.d;
    std::complex<ffloat> A_1_sigma=(u*u-1i*u)*.5*d_sigma*hlp.cosh_v;
    std::complex<ffloat> A_2_sigma=p.rho/p.sigma*A_2_rho+(2.+tau*hlp.xi)/(p.v_0*tau*hlp.xi*1i*u)*A_1_rho+p.sigma*tau*hlp.A_1/p.v_0*.5;
    std::complex<ffloat> A_sigma=A_1_sigma/hlp.A_2-hlp.A/hlp.A_2*A_2_sigma;
    
    std::complex<ffloat> tiuvmdivs=p.v_m*tau*1i*u/p.sigma;
    std::complex<ffloat> sp2=p.sigma*p.sigma;
    std::complex<ffloat> kvm2divsp2=2.*p.kappa*p.v_m/sp2;
    
    return{chf_val*(-hlp.A/p.v_0),
            chf_val*(2.*p.kappa/sp2*hlp.D+p.kappa*p.rho*tau*1i*u/p.sigma),
            chf_val*(-A_rho+kvm2divsp2/hlp.d*(d_rho-hlp.d/hlp.A_2*A_2_rho)+tiuvmdivs*p.kappa),
            chf_val*(-A_rho/(p.sigma*1i*u)+2.*p.v_m/sp2*hlp.D+kvm2divsp2/hlp.B*B_kappa+tiuvmdivs*p.rho),
            chf_val*(-A_sigma-2.*kvm2divsp2/p.sigma*hlp.D+kvm2divsp2/hlp.d*(d_rho-hlp.d/hlp.A_2*A_2_sigma)-tiuvmdivs/p.sigma*p.rho*p.kappa)};
}
ffloat HDistribution::first_order_moment()
{
    return -0.5 * p.v_m * T;
}
ffloat HDistribution::second_order_moment()
{
    auto const& [kappa, s2, k, r, v] = p;
    ffloat kappa2 = kappa * kappa, kappa3 = kappa2 * kappa;
    ffloat k2 = k * k;
    ffloat t = T;
    return s2 / (8 * kappa3) * (-k2 * std::exp(-2 * kappa * t) + 4 * k * std::exp(- kappa * t) * (k - 2 * kappa * r) + 2 * kappa * t * (4 * kappa2 + k2 - 4 * kappa * k * r) + k * (8 * kappa * r - 3 * k));
}
ffloat HDistribution::fourth_order_moment()
{
    auto const& [a, s2, k, r, v] = p;
    ffloat a2 = a * a, a3 = a2 * a, a4 = a3 * a;
    ffloat k2 = k * k, k3 = k2 * k, k4 = k3 * k;
    ffloat t = T, t2 = t * t;
    ffloat r2 = r * r;
    return (3 * k2 * s2) / (64 * std::pow(a, 7)) * (
        -3 * k4 * std::exp(-4 * a * t)
        -8 * k2 * std::exp(-3 * a * t) * (2 * a * k * t * (k - 2 * a * r) + 4 * a2 + k2 - 6 * a * k * r)
        -4 * std::exp(-2 * a * t) * (4 * a2 * k2 * t2 * std::pow(k - 2 * a * r, 2) + 2 * a * k * t * (k3 - 16 * a3 * r -12 * a * k2 * r + 4 * a2 * k * (3 + 4 * r2)) + 8 * a4 - 3 * k4 - 32 * a3 * k * r + 8 * a * k3 * r + 16 * a2 * k2 * r2)
        -8 * std::exp(-a * t) * (- 2 * a2 * k * t2 * std::pow(k - 2 * a * r, 3) - 8 * a * t * (k4 - 7 * a * k3 * r + 4 * a4 * r2 - 8 * a3 * k * r * (1 + r2) + a2 * k2 * (3 + 14 * r2)) - 9 * k4 + 70 * a * k3 * r + 32 * a3 * k * r * (4 + 3 * r2) - 16 * a4 * (1 + 4 * r2) - 4 * a2 * k2 * (9 + 40 * r2))
        + 4 * a* t * (5 * k4 - 40 * a * k3 * r - 32 * a3 * k * r * (3 + 2 * r2) + 16 * a4 * (1 + 4 * r2) + 24 * a2 * k2 * (1 + 4 * r2))
        -73 * k4 + 544 * a * k3 * r + 128 * a3 * k * r * (7 + 6 * r2) - 32 * a4 * (3 + 16 * r2) - 64 * a2 * k2 * (4 + 19 * r2));
}

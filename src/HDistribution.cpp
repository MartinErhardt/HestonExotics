#include"HDistribution.h"
#include"BSM.h"
#include <cmath>
#include<complex.h>
#include<cstdlib>
#include<iostream>
using namespace std::complex_literals;
ffloat HDistribution::int_error(const unsigned int trunc_m)const {
    ffloat exp2_m=std::exp2(trunc_m);
    return std::fabs(chf(exp2_m*M_PI)+chf(-exp2_m*M_PI))/(4*exp2_m*M_PI*M_PI*tau);
}
HDistribution::Helpers::Helpers(const HParams& p,const std::complex<ffloat> u,const ffloat tau){
        xi=p.kappa+p.sigma*p.rho*1i*u;
        std::complex<ffloat> fac=u*u-1i*u;
        d=std::sqrt(xi*xi+p.sigma*p.sigma*fac); //Why does this work?
        sinh_v=std::sinh(d*tau*.5); //numerically instable
        cosh_v=std::cosh(d*tau*.5);
        A_1=fac*sinh_v;
        A_2=d/p.v_0*cosh_v+xi/p.v_0*sinh_v;
        A=A_1/A_2;
        exp_kappa_tau=std::exp(p.kappa*tau*.5);
        B=d*exp_kappa_tau/(p.v_0*A_2);
        D=std::log((2.*d)/(d+xi+(d-xi)*std::exp(-d*tau)))+(p.kappa-d)*tau*.5;
        //if(d!=d||A_1!=A_1||A)
        //std::complex<ffloat> D2=std::log(d/p.v_0)+(p.kappa-d)*tau*.5-std::log((d+xi)/(2*p.v_0)+(d-xi)/(2*p.v_0)*std::exp(-d*tau));
        
}
std::complex<ffloat> HDistribution::chf(const std::complex<ffloat> u) const{
    HDistribution::helpers hlp(this->p,u,tau);
    return chf(u,hlp);
}
std::complex<ffloat> HDistribution::chf(const std::complex<ffloat> u,const helpers& hlp) const{
    std::complex<ffloat> res=std::exp((p.kappa*p.v_m*p.rho*tau*u*1i)/p.sigma-hlp.A+(2.*hlp.D*p.kappa*p.v_m)/(p.sigma*p.sigma)); //FIXME?
    //std::cout<<"1xi: "<<hlp.xi<<"\td: "<<hlp.d<<"\tA_1: "<<hlp.A_1<<"\tA_2: "<<hlp.A_2<<"\tA: "<<hlp.A<<"\tB: "<<hlp.B<<"\tD: "<<hlp.D<<"\tchf: "<<res<<'\n';
    //std::cout<<"res: "<<res<<'\n';
    //std::complex<double> char_u = exp(-v_bar * A_over_v0 + two_kappa_v_bar_over_sigma_squared * D - kappa_v_bar_rho_T_over_sigma * ui);
    return res;
}

std::vector<std::complex<ffloat>> HDistribution::chf_grad(const std::complex<ffloat> u) const{
    HDistribution::helpers hlp(this->p,u,tau);
    std::complex<ffloat> chf_val_arg=chf(u,hlp);
    return chf_grad(u,hlp,chf_val_arg);
}
std::vector<std::complex<ffloat>> HDistribution::chf_grad(const std::complex<ffloat> u,const helpers& hlp,const std::complex<ffloat>chf_val) const{
    std::vector<std::complex<ffloat>> chf_and_chf_grad = chf_chf_grad(u,hlp,chf_val);
    return {chf_and_chf_grad[1],chf_and_chf_grad[2],chf_and_chf_grad[3],chf_and_chf_grad[4],chf_and_chf_grad[5]};
}
std::vector<std::complex<ffloat>> HDistribution::chf_chf_grad(const std::complex<ffloat> u) const{
    HDistribution::helpers hlp(this->p,u,tau);
    std::complex<ffloat> chf_val_arg=chf(u,hlp);
    return chf_chf_grad(u,hlp,chf_val_arg);
}
std::vector<std::complex<ffloat>> HDistribution::chf_chf_grad(const std::complex<ffloat> u,const helpers& hlp,const std::complex<ffloat>chf_val) const{
    std::complex<ffloat> d_rho=hlp.xi*p.sigma*1i*u/hlp.d;
    std::complex<ffloat> A_2_rho=p.sigma*1i*u*(2.+hlp.xi*tau)/(2.*hlp.d*p.v_0)*(hlp.xi*hlp.cosh_v+hlp.d*hlp.sinh_v);
    std::complex<ffloat> B_rho=hlp.exp_kappa_tau/p.v_0*(d_rho/hlp.A_2-hlp.d/(hlp.A_2*hlp.A_2)*A_2_rho);
    std::complex<ffloat> A_1_rho=(1i*u*(u*u-1i*u)*tau*hlp.xi*p.sigma)/(2.*hlp.d)*hlp.cosh_v;
    std::complex<ffloat> A_rho=A_1_rho/hlp.A_2-hlp.A/hlp.A_2*A_2_rho;
    //std::complex<ffloat> A_kappa=-1i/(p.sigma*u)*A_rho;
    std::complex<ffloat> B_kappa=-1i/(p.sigma*u)*B_rho+hlp.B*tau*.5;
    std::complex<ffloat> d_sigma=(p.rho/p.sigma-1./hlp.xi)*d_rho+p.sigma*u*u/hlp.d;
    std::complex<ffloat> A_1_sigma=(u*u-1i*u)*.5*tau*d_sigma*hlp.cosh_v;
    std::complex<ffloat> A_2_sigma=p.rho/p.sigma*A_2_rho+(2.+tau*hlp.xi)/(p.v_0*tau*hlp.xi*1i*u)*A_1_rho+p.sigma*tau*hlp.A_1/p.v_0*.5;
    std::complex<ffloat> A_sigma=A_1_sigma/hlp.A_2-hlp.A/hlp.A_2*A_2_sigma;
    
    std::complex<ffloat> tiuvmdivs=p.v_m*tau*1i*u/p.sigma;
    std::complex<ffloat> sp2=p.sigma*p.sigma;
    std::complex<ffloat> kvm2divsp2=2.*p.kappa*p.v_m/sp2;
    
    std::complex<ffloat> h_v_0=(-hlp.A/p.v_0);
    std::complex<ffloat> h_v_m=(2.*p.kappa/sp2*hlp.D+p.kappa*p.rho*tau*1i*u/p.sigma);
    std::complex<ffloat> h_sigma=(-A_sigma-2.*kvm2divsp2/p.sigma*hlp.D+kvm2divsp2/hlp.d*(d_sigma-hlp.d/hlp.A_2*A_2_sigma)-tiuvmdivs/p.sigma*p.rho*p.kappa);
    std::complex<ffloat> h_kappa=(-A_rho/(p.sigma*1i*u)+2.*p.v_m/sp2*hlp.D+kvm2divsp2/hlp.B*B_kappa+tiuvmdivs*p.rho);
    std::complex<ffloat> h_rho=(-A_rho+kvm2divsp2/hlp.d*(d_rho-hlp.d/hlp.A_2*A_2_rho)+tiuvmdivs*p.kappa);
    
    return { chf_val,
            chf_val*h_v_0,
            chf_val*h_v_m,
            chf_val*h_rho,
            chf_val*h_kappa,
            chf_val*h_sigma
           };
}

ffloat HDistribution::first_order_moment() const
{
    return -.5 * p.v_m * tau;
}

ffloat HDistribution::second_order_moment() const
{
    auto const& [v,s2,r,kappa,k] = p;
    ffloat kappa2 = kappa * kappa, kappa3 = kappa2 * kappa;
    ffloat k2 = k * k;
    ffloat t = tau;
    return s2 / (8 * kappa3) * (-k2 * std::exp(-2 * kappa * t) + 4 * k * std::exp(- kappa * t) * (k - 2 * kappa * r) + 2 * kappa * t * (4 * kappa2 + k2 - 4 * kappa * k * r) + k * (8 * kappa * r - 3 * k));
}
ffloat HDistribution::fourth_order_moment() const
{
    auto const& [v,s2,r,a,k] = p;
    ffloat a2 = a * a, a3 = a2 * a, a4 = a3 * a;
    ffloat k2 = k * k, k3 = k2 * k, k4 = k3 * k;
    ffloat t = tau, t2 = t * t;
    ffloat r2 = r * r;
    return (3 * k2 * s2) / (64 * std::pow(a, 7)) * (
        -3 * k4 * std::exp(-4 * a * t)
        -8 * k2 * std::exp(-3 * a * t) * (2 * a * k * t * (k - 2 * a * r) + 4 * a2 + k2 - 6 * a * k * r)
        -4 * std::exp(-2 * a * t) * (4 * a2 * k2 * t2 * std::pow(k - 2 * a * r, 2) + 2 * a * k * t * (k3 - 16 * a3 * r -12 * a * k2 * r + 4 * a2 * k * (3 + 4 * r2)) + 8 * a4 - 3 * k4 - 32 * a3 * k * r + 8 * a * k3 * r + 16 * a2 * k2 * r2)
        -8 * std::exp(-a * t) * (- 2 * a2 * k * t2 * std::pow(k - 2 * a * r, 3) - 8 * a * t * (k4 - 7 * a * k3 * r + 4 * a4 * r2 - 8 * a3 * k * r * (1 + r2) + a2 * k2 * (3 + 14 * r2)) - 9 * k4 + 70 * a * k3 * r + 32 * a3 * k * r * (4 + 3 * r2) - 16 * a4 * (1 + 4 * r2) - 4 * a2 * k2 * (9 + 40 * r2))
        + 4 * a* t * (5 * k4 - 40 * a * k3 * r - 32 * a3 * k * r * (3 + 2 * r2) + 16 * a4 * (1 + 4 * r2) + 24 * a2 * k2 * (1 + 4 * r2))
        -73 * k4 + 544 * a * k3 * r + 128 * a3 * k * r * (7 + 6 * r2) - 32 * a4 * (3 + 16 * r2) - 64 * a2 * k2 * (4 + 19 * r2));
}

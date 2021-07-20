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
    //std::cout<<
    //"d_rho: "<<d_rho<<"\tA2_rho: "<<A_2_rho<<"\tB_rho: "<<B_rho<<"\tA1_rho: "<<A_1_rho<<"\tA_rho: "<<A_rho<<
    //"\td_sigma: "<<d_sigma<<"\tA_1_sigma: "<<A_1_sigma<<"\tA_2_sigma: "<<A_2_sigma<<"\tA_sigma: "<<A_sigma<<'\n';
    //std::cout<<"MY char_u: "<<chf_val<<"\th_v_0: "<<h_v_0<<"\th_v_bar: "<< h_v_m<<"\th_rho: "<<h_rho<<"\th_kappa: "<<h_kappa<<"\th_sigma: "<<h_sigma<<'\n';
    return { chf_val,
            chf_val*h_v_0,
            chf_val*h_v_m,
            chf_val*h_rho,
            chf_val*h_kappa,
            chf_val*h_sigma
           };
}
std::vector<std::complex<double>> GetCuiGradient(HParams p, std::complex<double> u, double T)
{
    auto const& [v0, v_bar, rho,kappa, sigma] = p;
    double sigma_times_rho=sigma*rho;
    double sigma_squared=sigma*sigma;
    double half_T=T*.5;
    double two_kappa_v_bar_over_sigma_squared=2*kappa*v_bar/(sigma*sigma);
    double kappa_v_bar_rho_T_over_sigma=kappa*v_bar*rho*T/sigma;
    
    std::complex<double> ui = 1i * u;
    std::complex<double> u_squared = u * u;

    std::complex<double> xi = kappa - sigma_times_rho * ui; // xi = kappa - sigma * rho * u * i
    std::complex<double> xi2 = xi * xi;
    std::complex<double> m = ui + u_squared; // m = u * i + u^2;
    std::complex<double> d = std::sqrt(xi2 + m * sigma_squared); // d = sqrt(pow(xi,2) + m*pow(sigma,2));

    // alp, calp, salp
    std::complex<double> alpha = d * half_T;
    std::complex<double> cosh_alpha = cosh(alpha);
    std::complex<double> sinh_alpha = sinh(alpha);
    std::complex<double> A2_times_v0 = d * cosh_alpha + xi * sinh_alpha;
    std::complex<double> A1 = m * sinh_alpha;
    std::complex<double> A_over_v0 = A1 / A2_times_v0;

    std::complex<double> D = std::log(d) + (kappa - d) * half_T - std::log((d + xi) * 0.5 + (d - xi) * 0.5 * exp(-d * T));
    
//        std::complex<double> g = exp(minus_kappa_v_bar_rho_T_over_sigma * ui);

    // F = S * e^((r - q) * T)
    // characteristic function: y1 = exp(i * log(F / S) * u) * exp(-A + 2 * kappa * b / pow(sigma, 2) * D) * g
    // But we only care about the second exponential, the rest depends only on market parameters and will be computed separately.
    std::complex<double> char_u =std::exp(-v0 * A_over_v0 + two_kappa_v_bar_over_sigma_squared * D - kappa_v_bar_rho_T_over_sigma * ui);

    // B = d * exp(kappa * T / 2) / (A2 * v0);
    double exp_kappa_times_half_T = exp(kappa * half_T); // exp(kappa * T / 2)
    std::complex<double> B = d * exp_kappa_times_half_T / A2_times_v0;
    //std::cout<<"2xi: "<<xi<<"\td: "<<d<<"\tA_1: "<<A1<<"\tA_2: "<<A2_times_v0/v0<<"\tA"<<A_over_v0*v0<<"\tB: "<<B<<"\tD: "<<D<<"\tchf: "<<char_u<<'\n';
    // g = exp(-kappa * b * rho * T * u1 * i / sigma);
    double kappa_v_bar_rho_T = kappa * v_bar * rho * T;  // TAG: PRECOMPUTE
    double minus_kappa_v_bar_rho_T_over_sigma = -kappa_v_bar_rho_T / sigma;  // TAG: PRECOMPUTE

    std::complex<double> H = xi * cosh_alpha + d * sinh_alpha;

    // lnB = log(B);
    std::complex<double> lnB = D;

    // partial b: y3 = y1*(2*kappa*lnB/pow(sigma,2)-kappa*rho*T*u1*i/sigma);
    double two_kappa_over_sigma_squared = two_kappa_v_bar_over_sigma_squared / v_bar;
    double minus_kappa_rho_T_over_sigma = minus_kappa_v_bar_rho_T_over_sigma / v_bar;

    std::complex<double> h_v_bar = two_kappa_over_sigma_squared * lnB + minus_kappa_rho_T_over_sigma * ui;

    // partial rho:
    double minus_kappa_v_bar_t_over_sigma = minus_kappa_v_bar_rho_T_over_sigma / rho; //-kappa * v_bar * T/sigma;

    std::complex<double> sigma_ui_over_d = sigma * ui / d;
    std::complex<double> pd_prho = -xi * sigma_ui_over_d;
    std::complex<double> pA1_prho = m * cosh_alpha * half_T * pd_prho;
    std::complex<double> pA2_prho = -sigma_ui_over_d * H * (1.0 + xi * half_T);
    std::complex<double> pA_prho = (pA1_prho - A_over_v0 * pA2_prho) / A2_times_v0;
    std::complex<double> pd_phrho_minus_pA2_prho_times_d_over_A2 = pd_prho - pA2_prho * d / A2_times_v0;
    std::complex<double> pB_prho = exp_kappa_times_half_T / A2_times_v0 * pd_phrho_minus_pA2_prho_times_d_over_A2;
    std::complex<double> h_rho = -v0 * pA_prho + two_kappa_v_bar_over_sigma_squared * pd_phrho_minus_pA2_prho_times_d_over_A2 / d + minus_kappa_v_bar_t_over_sigma * ui;

    // partial kappa:
    double v_bar_rho_T_over_sigma = v_bar * rho * T / sigma;
    double two_v_bar_over_sigma_squared = two_kappa_v_bar_over_sigma_squared / kappa; // 2 * v_bar / sigma_squared;

    std::complex<double> minus_one_over_sigma_ui = -1.0 / (sigma * ui);
    std::complex<double> pB_pa = minus_one_over_sigma_ui * pB_prho + B * half_T;
    std::complex<double> h_kappa = -v0 * pA_prho * minus_one_over_sigma_ui + two_v_bar_over_sigma_squared * lnB + kappa * two_v_bar_over_sigma_squared * pB_pa / B - v_bar_rho_T_over_sigma * ui;

    // partial sigma:
    double rho_over_sigma = rho / sigma;
    double four_kappa_v_bar_over_sigma_cubed = 4 * kappa * v_bar / pow(sigma, 3);
    double kappa_v_bar_rho_T_over_sigma_squared = kappa_v_bar_rho_T / sigma_squared;
    std::complex<double> pd_pc = (rho_over_sigma - 1.0 / xi) * pd_prho + sigma * u_squared / d;
    std::complex<double> pA1_pc = m * cosh_alpha * half_T * pd_pc;
    std::complex<double> pA2_pc = rho_over_sigma * pA2_prho - 1.0 / ui * (2.0 / (T * xi) + 1.0) * pA1_prho + sigma * half_T * A1;
    std::complex<double> pA_pc = pA1_pc / A2_times_v0 - A_over_v0 / A2_times_v0 * pA2_pc;
    std::complex<double> h_sigma = -v0 * pA_pc - four_kappa_v_bar_over_sigma_cubed * lnB + two_kappa_v_bar_over_sigma_squared / d * (pd_pc - d / A2_times_v0 * pA2_pc) + kappa_v_bar_rho_T_over_sigma_squared * ui;
    //std::cout<<
    //"d_rho: "<<pd_prho<<"\tA2_rho: "<<pA2_prho/v0<<"\tB_rho: "<<pB_prho<<"\tA1_rho: "<<pA1_prho<<"\tA_rho: "<<v0*pA_prho<<
    //"\td_sigma: "<<pd_pc<<"\tA_1_sigma: "<<pA1_pc<<"\tA_2_sigma: "<<pA2_pc/v0<<"\tA_sigma: "<<pA_pc*v0<<'\n';
    //std::cout<<"   char_u: "<<char_u<<"\th_v_0: "<<-A_over_v0<<"\th_v_bar: "<< h_v_bar<<"\th_rho: "<<h_rho<<"\th_kappa: "<<h_kappa<<"\th_sigma: "<<h_sigma<<'\n';
    return { -A_over_v0 * char_u, char_u * h_v_bar, char_u * h_rho, char_u * h_kappa,  char_u * h_sigma};
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
/*
	double m_kappa;     // mean reversion rate
	double m_v_bar;     // long term variance
	double m_sigma;     // variance of volatility
	double m_rho;       // correlation between spot and volatility
	double m_v0;        // initial variance
    ffloat v_0; // inital varince
    ffloat v_m; // long term variance
    ffloat rho; // correlation between spot and volatility
    ffloat kappa; // mean reversion rate
    ffloat sigma; // variance of volatility [kappa, s2, k, r, v]
*/
void distr_test(){

    /////// FX
    //double kappa = 0.5;           // |  mean reversion rate
    //double v_bar = 0.04;          // |  long term variance
    //double sigma = 1;          // |  variance of volatility
    //double rho = -0.9;            // |  correlation between spot and volatility
    //double v0 = 0.04;             // |  initial variance
    /////// IR
    //    double kappa = 0.3;           // |  mean reversion rate
    //    double v_bar = 0.04;          // |  long term variance
    //    double sigma = 0.9;          // |  variance of volatility
    //    double rho = -0.5;            // |  correlation between spot and volatility
    //    double v0 = 0.04;             // |  initial variance
    /////// EQ
        double kappa = 1;           // |  mean reversion rate
        double v_bar = 0.09;          // |  long term variance
        double sigma = 1;          // |  variance of volatility
        double rho = 0.04;            // |  correlation between spot and volatility
        double v0 = 0.09;             // |  initial variance

    /*
    for(int i=1;i<1000; i++){
        test.p.v_0=v0 * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.v_m=v_bar * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.rho=rho * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.kappa=kappa * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        test.p.sigma=sigma * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5);
        double tau=(((double) rand() / (RAND_MAX)));// *trading_days;
        double x=(((double) rand() / (RAND_MAX)) - 0.5)*100;
        std::complex<double> correct =chf2(test.p,x,tau);
        std::complex<double> wrong = test.chf(x,tau);
        std::cout<<"v_0: "<<test.p.v_0<<"\tv_m: "<<test.p.v_m<<"\trho: "<<test.p.rho<<"\tkappa: "<<test.p.kappa<<"\tsigma: "<<test.p.sigma<<"\ttau: "<<tau<<"\tx: "<<x<<"\tmy_chf: "<<wrong<<"\tEudald Romo chf: "<< correct<<"\tdiff: "<<std::fabs(correct-wrong)<<'\n';
    }*/
     for(int i=1;i<1000; i++){
         HDistribution test(
             {      v0 * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    v_bar * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    rho * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    kappa * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5),
                    sigma * (1.0 + (((double) rand() / (RAND_MAX)) - 0.5) / 5)
            }, (((double) rand() / (RAND_MAX))));
        double x=(((double) rand() / (RAND_MAX)) - 0.5)*1000;
        std::vector<std::complex<double>> correct =GetCuiGradient(test.p,-x,test.tau);
        std::vector<std::complex<double>> wrong = test.chf_grad(x);
        std::cout<<"v_0: "<<test.p.v_0<<"\tv_m: "<<test.p.v_m<<"\trho: "<<test.p.rho<<"\tkappa: "<<test.p.kappa<<"\tsigma: "<<test.p.sigma<<"\ttau: "<<test.tau<<"\tx: "<<x<<"\tdiff: "<<std::sqrt((correct[0]-wrong[0])*std::conj((correct[0]-wrong[0]))+
                    (correct[1]-wrong[1])*std::conj((correct[1]-wrong[1]))+
                    (correct[2]-wrong[2])*std::conj((correct[2]-wrong[2]))+
                    (correct[3]-wrong[3])*std::conj((correct[3]-wrong[3]))+
                    (correct[4]-wrong[4])*std::conj((correct[4]-wrong[4])))<<'\n';
    }
}

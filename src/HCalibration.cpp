/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HCalibration.h"
#include"BSM.h"
#include"SWIFT.h"
#include<iostream>
#include<levmar/levmar.h>
/*typedef struct {
    ffloat mu;
    ffloat theta;
    ffloat kappa;
    ffloat xi;
    ffloat rho;
} HestonParams;
*/
typedef std::numeric_limits< double > dbl;

typedef struct ED{
    const options_chain& opts;
    HDistribution* distr;
    std::shared_ptr<SWIFT> pricing_method=nullptr;
    ED(const options_chain& opts,HParams p, ffloat expi): opts(opts),distr(new HDistribution(p,expi)){}
    ~ED(){ //std::cout<<"ED deleted\n"; 
        //delete distr;
    } //moved when push_back 
} expiry_data;
typedef struct AS{
    ffloat S;
    std::vector<expiry_data>& exp_list;
    ~AS(){//std::cout<<"freed? "<<'\n';
        for(auto e :exp_list){
        //std::cout<<"freed? "<<e.distr<<'\n';
        delete (e.distr);
    }
    } //moved when push_back 
} adata_s;
#define NEWOLD_METHOD_RATIO 0.1
void update_adata(ffloat *p, adata_s * adata);
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata);
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata);

void update_adata(ffloat *p, adata_s * adata){
    const HParams new_params={p[0],p[1],p[2],p[3],p[4]};
    //std::shared_ptr<SWIFT> current=nullptr;
    for(auto exp_data : adata->exp_list){
        if(!(exp_data.distr->p==new_params)){
            exp_data.pricing_method->flush_cache();
            //std::cout<<"delete distr\n";
            exp_data.distr->p=new_params;
            //std::cout<<"new params at: "<<&(exp_data.distr->p)<<"\tv_0: "<<exp_data.distr->p.v_0<<"\tv_m: "<<exp_data.distr->p.v_m<<"\trho: "<<exp_data.distr->p.rho<<"\tkappa"<<exp_data.distr->p.kappa<<"\tsigma: "<<exp_data.distr->p.sigma<<'\n';
        }
        //std::cout<<"cur params at: "<<&(exp_data.distr->p)<<"\tv_0: "<<(float)exp_data.distr->p.v_0<<"\tv_m: "<<exp_data.distr->p.v_m<<"\trho: "<<exp_data.distr->p.rho<<"\tkappa"<<exp_data.distr->p.kappa<<"\tsigma: "<<exp_data.distr->p.sigma<<'\n';
        
        auto new_swift_parameters=SWIFT::get_parameters(*exp_data.distr,adata->S,exp_data.opts);
        if (new_swift_parameters->m>exp_data.pricing_method->my_params.m 
            //||new_swift_parameters->J<NEWOLD_METHOD_RATIO*exp_data.pricing_method->my_params.J
        ){
            //delete exp_data.pricing_method;
        //    if(current==nullptr||new_swift_parameters->m>current->my_params.m ||new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J)
        //        current=std::make_shared<SWIFT>(*new_swift_parameters);
            exp_data.pricing_method=std::make_shared<SWIFT>(*new_swift_parameters);//std::shared_ptr(current);
        }
    }
}
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata){
    std::cout<<"get jac\t\t v_0: "<<p[0]<<"\tv_m: "<<p[1]<<"\trho: "<<p[2]<<"\tkappa: "<<p[3]<<"\tsigma: "<<p[4]<<'\n';
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* jac2=jac;
    //std::cout<<"dimension m: "<<m<<"\tdimension n: "<<n_observations<<"\tbuffer start"<<jac2<<"\tbuffer end: "<<jac+n_observations*m<<"\tbuffer size: "<<jac+n_observations*m-jac2<<'\n';
    for(auto exp_data : my_adata->exp_list) exp_data.pricing_method->price_opts_grad(*exp_data.distr,my_adata->S,exp_data.opts, &jac2,jac+n_observations*m);
    if(jac2<jac+n_observations*m) throw std::runtime_error("Gradient buffer too large");
}
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata){
    std::cout.precision(dbl::max_digits10);
    std::cout<<"get prices\t v_0: "<<p[0]<<"\tv_m: "<<p[1]<<"\trho: "<<p[2]<<"\tkappa: "<<p[3]<<"\tsigma: "<<p[4]<<'\n';
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* x2=x;
    //std::cout<<"dimension m: "<<m<<"\tdimension n: "<<n_observations<<"\tbuffer start"<<x2<<"\tbuffer end: "<<x+n_observations<<"\tbuffer size: "<<x+n_observations-x2<<'\n';
    for(auto exp_data : my_adata->exp_list){
        //ffloat* x3=x2;
        exp_data.pricing_method->price_opts(*exp_data.distr,my_adata->S,exp_data.opts, &x2,x+n_observations);
        //std::cout<<"Wrote x bytes to buffer: "<<x3-x2<<'\n';
    }
    if(x2<x+n_observations) throw std::runtime_error("Pricing buffer too large");
    //for(ffloat* x2=x;x2<x+n_observations;x2++) std::cout<<*x2<<'\n';
}
std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data){
    adata_s adata={S,*(new std::vector<expiry_data>())};
    //HDistribution current_distribution;
    //std::shared_ptr<SWIFT> current;
    ffloat avg_iv=avg_imp_vol(S, market_data);
    // We start in a Black-Scholes Model see p.106 lecture notes
    ffloat yearly_avg_iv = pow(1.+avg_iv*avg_iv,trading_days)-1.;
    ffloat p[5];
    p[0]=yearly_avg_iv;
    p[1]=yearly_avg_iv;
    p[2]=.04;
    p[3]=1.;
    p[4]=1.;
    unsigned int n_observations_cur=0;
    for(std::list<options_chain>::const_iterator opts = market_data.begin(); opts != market_data.end(); opts++){
        if(!opts->days_to_expiry) continue;
        if(opts->options.size()){
            n_observations_cur+=opts->options.size();
            expiry_data *expi_cur=new expiry_data(*opts,{p[0],p[1],p[2],p[3],p[4]},opts->time_to_expiry);
            auto new_swift_parameters=SWIFT::get_parameters(*expi_cur->distr,S,*opts);
            //if (current==nullptr|| new_swift_parameters->m>current->my_params.m || new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J){
                //std::cout<<"is here the free1?\n";
                //current=std::make_shared<SWIFT>(*new_swift_parameters);
                //std::cout<<"is here the free?\n";
            //}
            expi_cur->pricing_method=std::make_shared<SWIFT>(*new_swift_parameters);//std::shared_ptr(current);
            adata.exp_list.push_back(*expi_cur);
        }
    }
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*n_observations_cur);
    ffloat * x2=x;
    for(std::list<options_chain>::const_iterator opts = market_data.begin(); opts != market_data.end(); opts++){
        if(!opts->days_to_expiry) continue;
        for(auto e :opts->options) *(x2++)=e.price; //test
    }
    std::cout<<"setup completed\t# observations: "<<n_observations_cur<<'\n';
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    int retval;
    if((retval=dlevmar_der(get_prices_for_levmar, get_jacobian_for_levmar, p, x, 5, n_observations_cur, 100, opts, info, NULL, NULL, (void*) &adata))<0) throw std::runtime_error("levmar failed!");
    auto to_calib=std::unique_ptr<HParams>(new HParams({p[0],p[1],p[2],p[3],p[4]}));
    std::cout<<"# iter: "<<retval<<"\tv_0: "<<to_calib->v_0<<"\tv_m: "<<to_calib->v_m<<"\trho: "<<to_calib->rho<<"\tkappa: "<<to_calib->kappa<<"\tsigma: "<<to_calib->sigma<<'\n';
    //std::cout<<"x: "<<x<<"\tx2: "<<x2<<'\n';
    
    //std::cout<<"double free here?\n";
    //delete &adata;
    return to_calib;
}

void pricing_test(){
    adata_s adata={1.,*(new std::vector<expiry_data>())};
    std::vector<ffloat> prices={
        0.079676812094469612,
        0.042586263756033402,
        0.023642876097251266,
        0.0019447635553313004,
        0.00022457675334788794,
        0.15394308829999556,
        0.062423699148410512,
        0.03716156027004304,
        0.0075329906749080954,
        0.0021671474663124877,
        0.20438572427747337,
        0.081424876885654141,
        0.047007477654992851,
        0.013394276263081869,
        0.0051665880568379507,
        0.24157958905122354,
        0.09899570130848076,
        0.055258126136188711,
        0.017997398670040403,
        0.0084472770368514277,
        0.27251527545827303,
        0.11480331435058673,
        0.062447187930262341,
        0.023683757123971076,
        0.011714319060476305,
        0.29968779707061638,
        0.12892214079993489,
        0.068994812796204188,
        0.02829225908895841,
        0.01488391596622375,
        0.3569911166621863,
        0.16022472995208151,
        0.086978384558870039,
        0.041295845512103642,
        0.024161759165286907,
        0.41608092405221836,
        0.18693437912935496,
        0.10163212168506526,
        0.052461532814837883,
        0.032268432242168528};
    std::vector<double> expiries = {
    0.119047619047619,
    0.238095238095238,
    0.357142857142857,
    0.476190476190476,
    0.595238095238095,
    0.714285714285714,
    1.07142857142857,
    1.42857142857143};
    std::vector<std::vector<double>> K_over_S = {
        {0.9371, 0.9956, 1.0427, 1.2287, 1.3939},
        {0.8603, 0.9868, 1.0463, 1.2399, 1.4102},
        {0.8112, 0.9728, 1.0499, 1.2485, 1.4291},
        {0.7760, 0.9588, 1.0530, 1.2659, 1.4456},
        {0.7470, 0.9464, 1.0562, 1.2646, 1.4603},
        {0.7216, 0.9358, 1.0593, 1.2715, 1.4736},
        {0.6699, 0.9175, 1.0663, 1.2859, 1.5005},
        {0.6137, 0.9025, 1.0766, 1.3046, 1.5328}
    };
    std::vector<swift_parameters> params={{5, 32, 5.6568542494923806, -1.8343994686679572, 1.5720210785285174, -58, 50, 128},
                                            {5, 32, 5.6568542494923806, -2.6596460616572828, 2.4759124462331239, -85, 79, 256},
                                            {5, 32, 5.6568542494923806, -3.3440060583236764, 3.2104875432489965, -107, 102, 256},
                                            {5, 32, 5.6568542494923806, -3.9452944230550919, 3.8494203406057474, -126, 123, 256},
                                            {5, 32, 5.6568542494923806, -4.4898573506150337, 4.4267150742095884, -143, 141, 256},
                                            {5, 32, 5.6568542494923806, -4.9920925407415853, 4.9592398930278314, -159, 158, 256},
                                            {5, 32, 5.6568542494923806, -6.3157233304280682, 6.353408918854841, -202, 203, 512},
                                            {5, 32, 5.6568542494923806, -7.4606624599663149, 7.5789582584617987, -238, 242, 512}};
    double kappa = 1;           // |  mean reversion rate
    double v_bar = 0.09;          // |  long term variance
    double sigma = 1;          // |  variance of volatility
    double rho = 0.04;            // |  correlation between spot and volatility
    double v0 = 0.09;
    double p[5];
    p[0]=v0;p[1]=v_bar,p[2]=rho;p[3]=kappa;p[4]=sigma;
    for (unsigned int index = 0; index < expiries.size(); ++index)
    {
        std::vector<ffloat> cur=K_over_S[index];
        options_chain * opt_chain=new options_chain(static_cast<unsigned int>(expiries[index])/trading_days,expiries[index]);
        for(auto c: cur){
            option * new_opt =new option();
            new_opt->volume=1;
            new_opt->strike=c;
            new_opt->price=1.;
            opt_chain->options.push_back(*new_opt);
        }
        opt_chain->min_strike=cur[0];
        opt_chain->max_strike=cur[cur.size()-1];
        expiry_data *expi_cur=new expiry_data(*opt_chain,{v0,v_bar,rho,kappa,sigma},expiries[index]);
        expi_cur->pricing_method=std::make_shared<SWIFT>(params[index]);//std::shared_ptr(current);
        adata.exp_list.push_back(*expi_cur);
    }
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*40);
    get_prices_for_levmar(&p[0], x, 5, 40, (void *) &adata);
    for(int i=0;i<40;i++){
        std::cout<<"x: "<<x[i]<<"\tp: "<<prices[i]<<"\tdiff: "<<std::fabs(x[i]-prices[i])<<'\n';
    }
}

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
typedef struct ED{
    const options_chain& opts;
    HDistribution* distr;
    std::shared_ptr<SWIFT> pricing_method=nullptr;
    ED(const options_chain& opts,HParams p, ffloat expi): opts(opts),distr(new HDistribution(p,expi)){}
    //~ED(){std::cout<<"ED deleted\n"; delete distr;} moved when push_back 
} expiry_data;
typedef struct AS{
    ffloat S;
    std::vector<expiry_data>& exp_list;
    ~AS(){for(auto e :exp_list) delete e.distr;}
} adata_s;
#define NEWOLD_METHOD_RATIO 0.1
void update_adata(ffloat *p, adata_s * adata);
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata);
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata);

void update_adata(ffloat *p, adata_s * adata){
    const HParams new_params={p[0],p[1],p[2],p[3],p[4]};
    std::shared_ptr<SWIFT> current=nullptr;
    for(auto exp_data : adata->exp_list){
        if(!(exp_data.distr->p==new_params)){
            exp_data.pricing_method->flush_cache();
            delete exp_data.distr;
            exp_data.distr=new HDistribution(new_params, static_cast<ffloat>(exp_data.opts.days_to_expiry)/trading_days);
        }
        auto new_swift_parameters=SWIFT::get_parameters(*exp_data.distr,adata->S,exp_data.opts);
        if (new_swift_parameters->m>exp_data.pricing_method->my_params.m ||new_swift_parameters->J<NEWOLD_METHOD_RATIO*exp_data.pricing_method->my_params.J){
            //delete exp_data.pricing_method;
            if(current==nullptr||new_swift_parameters->m>current->my_params.m ||new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J)
                current=std::make_shared<SWIFT>(*new_swift_parameters);
            exp_data.pricing_method=std::shared_ptr(current);
        }
    }
}
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata){
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* jac2=jac;
    std::cout<<"dimension m: "<<m<<"\tdimension n: "<<n_observations<<'\n';
    for(auto exp_data : my_adata->exp_list){
        ffloat* jac3=jac2;
        exp_data.pricing_method->price_opts_grad(*exp_data.distr,my_adata->S,exp_data.opts, &jac2,jac+n_observations*m);
        std::cout<<"Wrote x bytes to buffer: "<<jac3-jac2<<"\ttotal bytes: "<<jac2-jac<<'\n';
    }
    if(jac2<jac+n_observations*m) throw std::runtime_error("Gradient buffer too large");
}
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata){
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* x2=x;
    std::cout<<"dimension m: "<<m<<"\tdimension n: "<<n_observations<<"\tbuffer start"<<x2<<"\tbuffer end: "<<x+n_observations*m<<"\tbuffer size: "<<x+n_observations*m-x2<<'\n';
    for(auto exp_data : my_adata->exp_list){
        ffloat* x3=x2;
        exp_data.pricing_method->price_opts(*exp_data.distr,my_adata->S,exp_data.opts, &x2,x+n_observations*m);
        std::cout<<"Wrote x bytes to buffer: "<<x3-x2<<'\n';
    }
    if(x2<x+n_observations*m) throw std::runtime_error("Pricing buffer too large");
}
std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data){
    adata_s adata={S,*(new std::vector<expiry_data>())};
    //HDistribution current_distribution;
    std::shared_ptr<SWIFT> current;
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
            expiry_data *expi_cur=new expiry_data(*opts,{p[0],p[1],p[2],p[3],p[4]},static_cast<ffloat>(opts->days_to_expiry)/trading_days);
            auto new_swift_parameters=SWIFT::get_parameters(*expi_cur->distr,S,*opts);
            if (current==nullptr|| new_swift_parameters->m>current->my_params.m || new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J){
                std::cout<<"is here the free1?\n";
                current=std::make_shared<SWIFT>(*new_swift_parameters);
                std::cout<<"is here the free?\n";
            }
            expi_cur->pricing_method=std::shared_ptr(current);
            adata.exp_list.push_back(*expi_cur);
        }
    }
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*n_observations_cur);
    for(std::list<options_chain>::const_iterator opts = market_data.begin(); opts != market_data.end(); opts++){
        if(!opts->days_to_expiry) continue;
        for(auto e :opts->options) *(x++)=e.price; //test
    }
    std::cout<<"setup completed n_observations: "<<n_observations_cur<<'\n';
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    dlevmar_der(get_prices_for_levmar, get_jacobian_for_levmar, p, x, 5, n_observations_cur, 100, opts, info, NULL, NULL, (void*) &adata);
    auto to_calib=std::unique_ptr<HParams>(new HParams({p[0],p[1],p[2],p[3],p[4]}));
    delete &adata.exp_list;
    return to_calib;
}

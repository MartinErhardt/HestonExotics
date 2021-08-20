/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HCalibration.h"
#include"BSM.h"
#include<iostream>
#include<levmar/levmar.h>
typedef std::numeric_limits< double > dbl;
#define NEWOLD_METHOD_RATIO 0.1
void update_adata(ffloat *p, adata_s * adata);
void update_adata(ffloat *p, adata_s * adata){
    const HParams new_params={p[0],p[1],p[2],p[3],p[4]};
    //std::shared_ptr<SWIFT> current=nullptr; delete
    for(auto &exp_data : adata->exp_list){
        if(!(exp_data.distr->p==new_params)){
            exp_data.pricing_method->flush_cache();
            //std::cout<<"delete distr\n";
            exp_data.distr->p=new_params;//}
            //std::cout<<"new params at: "<<&(exp_data.distr->p)<<"\tv_0: "<<exp_data.distr->p.v_0<<"\tv_m: "<<exp_data.distr->p.v_m<<"\trho: "<<exp_data.distr->p.rho<<"\tkappa"<<exp_data.distr->p.kappa<<"\tsigma: "<<exp_data.distr->p.sigma<<'\n';
        //std::cout<<"cur params at: "<<&(exp_data.distr->p)<<"\tv_0: "<<(float)exp_data.distr->p.v_0<<"\tv_m: "<<exp_data.distr->p.v_m<<"\trho: "<<exp_data.distr->p.rho<<"\tkappa"<<exp_data.distr->p.kappa<<"\tsigma: "<<exp_data.distr->p.sigma<<'\n';
        //std::cout<<"new SWIFT\n";
        std::unique_ptr<swift_parameters> new_swift_parameters=SWIFT::get_parameters(*exp_data.distr,adata->S,exp_data.opts);
        //if (new_swift_parameters->m>exp_data.pricing_method->my_params.m 
            //||new_swift_parameters->J<NEWOLD_METHOD_RATIO*exp_data.pricing_method->my_params.J
        //){
            //delete exp_data.pricing_method;
        //    if(current==nullptr||new_swift_parameters->m>current->my_params.m ||new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J)
        //        current=std::make_shared<SWIFT>(*new_swift_parameters);
            //std::cout<<"old m: "<<new_swift_parameters->m<<"\tnew m: "<<exp_data.pricing_method->my_params.m<<'\n'; 
            delete exp_data.pricing_method;
            exp_data.pricing_method=new SWIFT(*new_swift_parameters);//std::shared_ptr(current);
        }
    }
}
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata){
    std::cout<<"get gradient\tv_0: "<<p[0]<<"\tv_m: "<<p[1]<<"\trho: "<<p[2]<<"\tkappa: "<<p[3]<<"\tsigma: "<<p[4]<<'\n';
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* jac2=jac;
    for(auto &exp_data : my_adata->exp_list) exp_data.pricing_method->price_opts_grad(*exp_data.distr,my_adata->S,exp_data.opts, &jac2,jac+n_observations*m);
    if(jac2<jac+n_observations*m) throw std::runtime_error("Gradient buffer too large");
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata){
#pragma GCC diagnostic pop
    std::cout<<"get prices\tv_0: "<<p[0]<<"\tv_m: "<<p[1]<<"\trho: "<<p[2]<<"\tkappa: "<<p[3]<<"\tsigma: "<<p[4]<<'\n';
    adata_s * my_adata=static_cast<adata_s*>(adata);
    update_adata(p,my_adata);
    ffloat* x2=x;
    for(auto &exp_data : my_adata->exp_list){ exp_data.pricing_method->price_opts(*exp_data.distr,my_adata->S,exp_data.opts, &x2,x+n_observations);}
    if(x2<x+n_observations) throw std::runtime_error("Pricing buffer too large");
    //std::cout<<*my_adata;
    //unsigned int underpriced=0;
    //if(my_adata->real_prices) for(x2=x;x2<x+n_observations;x2++) if(*x2<my_adata->real_prices[x2-x]) underpriced++;
    //std::cout<<"share of underpriced: "<<static_cast<ffloat>(underpriced)/n_observations<<'\n';
}

std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data){
    adata_s adata={S,nullptr,*(new std::list<expiry_data>())};
    //HDistribution current_distribution;
    //std::shared_ptr<SWIFT> current;
    ffloat avg_iv=avg_imp_vol(S, market_data);
    // We start in a Black-Scholes Model see p.106 lecture notes
    ffloat yearly_avg_iv = avg_iv*avg_iv;
    //std::cout<<"yearly_avg_iv: "<<yearly_avg_iv<<'\n';
    ffloat p[5];
    p[0]=yearly_avg_iv;
    p[1]=yearly_avg_iv;
    p[2]=-.2;
    p[3]=1.;
    p[4]=1.;
    unsigned int n_observations_cur=0;
    std::cout<<"start levmar setup\n";
    for(const auto &opts: market_data){
        //if(opts->time_to_expiry<=EXP_LB) continue;
        if(opts.options->size()>0){
            n_observations_cur+=opts.options->size();
            HDistribution *new_distr=new HDistribution({p[0],p[1],p[2],p[3],p[4]},opts.time_to_expiry,0.0005);
            std::unique_ptr<swift_parameters> new_swift_parameters=SWIFT::get_parameters(*new_distr,S,opts);
            //if (current==nullptr|| new_swift_parameters->m>current->my_params.m || new_swift_parameters->J<NEWOLD_METHOD_RATIO*current->my_params.J){
                //std::cout<<"is here the free1?\n";
                //current=std::make_shared<SWIFT>(*new_swift_parameters);
                //std::cout<<"is here the free?\n";
            //}
            SWIFT* pricing_method=new SWIFT(*new_swift_parameters);//std::shared_ptr(current);
            adata.exp_list.emplace_back(opts,new_distr,pricing_method);
        }
    }
    //std::cout<<"allocate x\n";
    ffloat * x=(ffloat*) malloc(sizeof(ffloat)*(n_observations_cur+1));
    ffloat * x2=x;
    for(std::list<options_chain>::const_iterator opts = market_data.begin(); opts != market_data.end(); opts++) for(auto e :*opts->options) *(x2++)=e.price;
    adata.real_prices=x;
    std::cout<<"setup completed\t# observations: "<<n_observations_cur<<'\n';
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    std::vector<std::string> msg={{"1 -stopped by small gradient J^T e"},
                                {"2 - stopped by small Dp"},
                                {"3 - stopped by itmax"},
                                {"4 - singular matrix. Restart from current p with increased \\mu"},
                                {"5 - no further error reduction is possible. Restart with increased mu"},
                                {"6 - stopped by small ||e||_2"},
                                {"7 - stopped by invalid (i.e. NaN or Inf) \"func\" values; a user error"}};
    int iter=dlevmar_der(get_prices_for_levmar, get_jacobian_for_levmar, p, x, 5, n_observations_cur, 100, opts, info, NULL, NULL, (void*) &adata);
    if(iter<0) throw std::runtime_error("levmar failed!");
    auto to_calib=std::unique_ptr<HParams>(new HParams({p[0],p[1],p[2],p[3],p[4]}));
    std::cout<<"# iter: "<<iter<<"\tv_0: "<<to_calib->v_0<<"\tv_m: "<<to_calib->v_m<<"\trho: "<<to_calib->rho<<"\tkappa: "<<to_calib->kappa<<"\tsigma: "<<to_calib->sigma<<"\tinital e: "<<info[0]<<"\te: "<<info[1]<<"\treason: "<<msg[info[6]-1]<<'\n';
    if(info[6]!=6.&&info[6]!=2) throw std::runtime_error("levmar failed! ");
    return to_calib;
}

std::ostream& operator<<(std::ostream& out, adata_s const& as){
    //ffloat * cur_price=as.real_prices;
    for(auto& ed:as.exp_list){
        for(option& o: *ed.opts.options){
            out<<"S: "<<as.S//<<"\temp price"<<*(cur_price++)
            <<"\task: "<<o.price<<"\tbid: "<<o.bid<<"\tstrike: "<<o.strike<<"\tvolume: "<<o.volume<<"\tdays left: "<<ed.distr->tau<<std::endl;
        }
        std::cout<<ed.pricing_method->my_params;
    }
    return out;
}

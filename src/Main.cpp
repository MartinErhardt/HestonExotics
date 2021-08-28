/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include <iostream>
#ifndef WIN32
#include <unistd.h>
#endif
#include <curl/curl.h>
#include<list>
#include "UnitTest.h"
#include "WebAPI.h"
#include"BSM.h"
#include"Main.h"
#include"HCalibration.h"
#include"HSimulation.h"
#include"AsianContract.h"
#include<map>
const std::string token="RVjzAiRnplMr78OblRHnVOvmb2SA";
const std::map<std::string, program_option> arg_eval {
        { "-c", CALIBRATE },
        { "--calibrate", CALIBRATE},
        { "-p", PRICE },
        {"--price", PRICE},
        {"-d", DOWNLOAD},
        {"--download",DOWNLOAD},
        {"-t", TEST},
        {"--test",TEST},
        {"-h",HELP},
        {"--help", HELP}
};
#define DEFAULT_VOL_TYPE "open_interest"
#define DEFAULT_VOL 0
typedef struct UD{
    ffloat S;
    std::list<options_chain>* all_chains;
    ~UD(){delete all_chains;}
}underlying_data;
underlying_data get_ddata(int argc, char* argv[], int* cur_arg,char**vol,int* vol_n){
    if(*cur_arg==argc-3||*cur_arg==argc-1) throw SynopsisError();
    std::string underlying=std::string(argv[*cur_arg+1]);
    if(*cur_arg==argc-2||(*cur_arg<argc-2&&argv[*cur_arg+2][0]=='-')){
        *vol=(char*)DEFAULT_VOL_TYPE;
        *vol_n=DEFAULT_VOL;
        *cur_arg+=2;
    } else try{
        *vol=argv[*cur_arg+2];
        std::cout<<argv[*cur_arg+3]<<std::endl;
        *vol_n=std::stoi(argv[*cur_arg+3]);
        *cur_arg+=4;
    }catch(std::out_of_range const&){ throw SynopsisError();}
    auto Getter=std::make_unique<WebInterface::WebAPI>(token);
    return {Getter->get_stock_quote(underlying), Getter->get_all_option_chains(underlying,*vol,*vol_n)};
}
unsigned int length(const std::list<options_chain>& all_chains){
    unsigned int ret_val=0;
    for(const options_chain& opt_chain: all_chains) ret_val+=opt_chain.options->size();
    return ret_val;
}
int main(int argc, char *argv[]) {
    int cur_arg=1;
    char* vol;
    int vol_n;
    while(cur_arg<argc){
        auto iter=arg_eval.find(std::string(argv[cur_arg]));
        if(iter==arg_eval.end()) throw SynopsisError();
        switch(iter->second){
            case(CALIBRATE):{
                underlying_data ddata=get_ddata(argc,argv,&cur_arg,&vol,&vol_n);
                std::cout<<"Options data downloaded and parsed\n";
                calibrate(ddata.S,*ddata.all_chains);
                break;
            }
            case(PRICE):{
                if((cur_arg==argc-4||cur_arg==argc-6)&&std::string(argv[cur_arg+2])=="all"&&std::string(argv[cur_arg+1])=="asian"){
                    cur_arg+=2;
                    underlying_data ddata=get_ddata(argc,argv,&cur_arg,&vol,&vol_n);
                    HSimulation::PricingTool<ffloat,AAsianCallNonAdaptive>my_pricing_tool(1);
                    std::vector<ffloat>& results=*my_pricing_tool.price(calibrate(ddata.S,*ddata.all_chains),ddata.S,*ddata.all_chains,1,length(*ddata.all_chains),1e+2);
                    unsigned int i=0;
                    for(const options_chain& opt_chain: *ddata.all_chains) for(const option& opt: *(opt_chain.options))
                        std::cout<<"S: "<<ddata.S<<"\tstrike: "<<opt.strike<<"\tbid: "<<opt.bid<<"\task: "<<opt.price<<"\tasian-option-price: "<<results[i++]<<"\tvolume: "<<opt.volume<<"\timp vol: "<<imp_vol(ddata.S,opt,opt_chain.time_to_expiry)<<"\tlb: "<<ddata.S-std::exp(-yearly_risk_free*opt_chain.time_to_expiry)*opt.strike<<"\texpiry time: "<<opt_chain.time_to_expiry*trading_days<<'\n';
                } else throw SynopsisError();
                break;
            }
            case(DOWNLOAD):{
                underlying_data ddata=get_ddata(argc,argv,&cur_arg,&vol,&vol_n);
                for(const options_chain& opt_chain: *ddata.all_chains) for(const option& opt: *(opt_chain.options))
                    std::cout<<"S: "<<ddata.S<<"\tstrike: "<<opt.strike<<"\tbid: "<<opt.bid<<"\task: "<<opt.price<<"\tvolume: "<<opt.volume<<"\timp vol: "<<imp_vol(ddata.S,opt,opt_chain.time_to_expiry)<<"\tlb: "<<ddata.S-std::exp(-yearly_risk_free*opt_chain.time_to_expiry)*opt.strike<<"\texpiry time: "<<opt_chain.time_to_expiry*trading_days<<'\n';
                break;
            }
            case(TEST):{
                auto iter2=test_eval.find(argv[cur_arg+1]);
                if(cur_arg==argc-1||iter2==test_eval.end()) throw SynopsisError();
                iter2->second();
                cur_arg+=2;
                break;
            }
            case(HELP):{
                std::cout<<help_info<<std::endl;
                std::cout<<"Author: Martin Erhardt (merhardt@rhrk.uni-kl.de)"<<std::endl;
                cur_arg++;
                break;
            }
        }
    }
    return 0;
}

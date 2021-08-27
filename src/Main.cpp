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
int main(int argc, char *argv[]) {
    int cur_arg=1;
    char* vol;
    int vol_n;
    while(cur_arg<argc-1){
        auto iter=arg_eval.find(std::string(argv[cur_arg]));
        if(iter==arg_eval.end()) throw SynapsisError();
        switch(iter->second){
            case(CALIBRATE):{
                if(cur_arg==argc-3||cur_arg==argc-1) throw SynapsisError();
                std::string underlying=std::string(argv[cur_arg+1]);
                if(cur_arg==argc-2||(cur_arg<argc-2&&argv[cur_arg+2][0]=='-')){
                    vol=(char*)DEFAULT_VOL_TYPE;
                    vol_n=DEFAULT_VOL;
                    cur_arg+=2;
                } else try{
                    vol=argv[cur_arg+2];
                    vol_n=std::stoi(argv[cur_arg+3]);
                    cur_arg+=4;
                }catch(std::out_of_range const&){ throw SynapsisError();}
                auto Getter=std::make_unique<WebInterface::WebAPI>(token);
                ffloat S=Getter->get_stock_quote(underlying);
                std::list<options_chain>* all_chains=Getter->get_all_option_chains(underlying,vol,vol_n);
                std::cout<<"Options data downloaded and parsed\n";
                calibrate(S,*all_chains);
                delete all_chains;
                break;
            }
            case(PRICE):{
                break;
            }
            case(DOWNLOAD):{
                if(cur_arg==argc-3||cur_arg==argc-1) throw SynapsisError();
                std::string underlying=std::string(argv[cur_arg+1]);
                if(cur_arg==argc-2||(cur_arg<argc-2&&argv[cur_arg+2][0]=='-')){
                    vol=(char*)DEFAULT_VOL_TYPE;
                    vol_n=DEFAULT_VOL;
                    cur_arg+=2;
                } else try{
                    vol=argv[cur_arg+2];
                    vol_n=std::stoi(argv[cur_arg+3]);
                    cur_arg+=4;
                }catch(std::out_of_range const&){ throw SynapsisError();}
                auto Getter=std::make_unique<WebInterface::WebAPI>(token);
                ffloat S=Getter->get_stock_quote(underlying);
                std::list<options_chain>* all_chains=Getter->get_all_option_chains(underlying,vol,vol_n);
                std::cout << std::fixed << std::setprecision(5) << std::setfill('0');
                for(const options_chain& opt_chain: *all_chains) for(const option& opt: *(opt_chain.options))
                std::cout<<"S: "<<S<<"\tstrike: "<<opt.strike<<"\tbid: "<<opt.bid<<"\task: "<<opt.price<<"\tvolume: "<<opt.volume<<"\timp vol: "<<imp_vol(S,opt,opt_chain.time_to_expiry)<<"\tlb: "<<S-std::exp(-yearly_risk_free*opt_chain.time_to_expiry)*opt.strike<<"\texpiry time: "<<opt_chain.time_to_expiry*trading_days<<'\n';
                delete all_chains;
                cur_arg+=2;
                break;
            }
            case(TEST):{
                auto iter2=test_eval.find(argv[cur_arg+1]);
                if(cur_arg==argc-1||iter2==test_eval.end()) throw SynapsisError();
                iter2->second();
                cur_arg+=2;
                break;
            }
            case(HELP):{
                cur_arg++;
                break;
            }
        }
    }
    return 0;
}

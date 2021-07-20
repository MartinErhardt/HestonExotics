/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include <iostream>
#ifndef WIN32
#include <unistd.h>
#endif
#include <curl/curl.h>
#include<list>
#include "HCalibration.h"
#include "HDistribution.h"
#include "WebAPI.h"
#define MY_TOKEN "RVjzAiRnplMr78OblRHnVOvmb2SA"
//ffloat call_price(const ffloat S, const ffloat K, const ffloat r, const ffloat sigma, const double T);
int main(int argc, char *argv[]) {
    std::string tmp_token(MY_TOKEN); // initialized on the stack s.t. reference is not temporary and comprimised.
    auto Getter=std::make_unique<WebInterface::WebAPI>(tmp_token);
    if (argc >=2 && std::string(*(argv+1)) == "quote"){
        ffloat S=Getter->get_stock_quote(std::string(*(argv+2)));
        //std::cout<<std::string(*(argv+2))<<"price: \t"<<S<<'\n';
        std::unique_ptr<std::list<options_chain>> all_chains=Getter->get_all_option_chains(std::string(*(argv+2)));
        std::cout<<"Options data downloaded and parsed\n";
        calibrate(S,*all_chains);
    }else if (argc ==3 && std::string(*(argv+1)) == "test" && std::string(*(argv+2)) == "distr"){
        distr_test();
    }else if (argc ==3 && std::string(*(argv+1)) == "test" && std::string(*(argv+2)) == "pricing"){
        pricing_test();
    }/*
    std::string input;
    std::cout<<"Enter S: ";
    std::getline (std::cin,input);
    ffloat S=std::stod(input);
    std::cout<<"Enter K: ";
    std::getline (std::cin,input);
    ffloat K=std::stod(input);
    std::cout<<"Enter r: ";
    std::getline (std::cin,input);
    ffloat r=std::stod(input);
    std::cout<<"Enter sigma: ";
    std::getline (std::cin,input);
    ffloat sigma=std::stod(input);
    std::cout<<"Enter time: ";
    std::getline (std::cin,input);
    double T=std::stod(input);
    std::cout<<"call price: "<< call_price(S,K,r,sigma,T)<<'\n';*/
    return 0;
}

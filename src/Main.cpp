/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include <iostream>
#ifndef WIN32
#include <unistd.h>
#endif
#include <curl/curl.h>

#include "SDE.h"
#include "WebAPI.h"
//#include "TWS.h"
#define MY_TOKEN "RVjzAiRnplMr78OblRHnVOvmb2SA"

int main(int argc, char *argv[]) {
    std::string tmp_token(MY_TOKEN); // initialized on the stack s.t. reference is not temporary and comprimised.
    auto Getter=std::make_unique<WebInterface::WebAPI>(tmp_token); 
    if (argc >=2 && std::string(*(argv+1)) == "quote"){
        ffloat S=Getter->get_stock_quote(std::string(*(argv+2)));
        //std::cout<<std::string(*(argv+2))<<"price: \t"<<S<<'\n';
        auto options_list=Getter->get_all_option_chains(std::string(*(argv+2)));
    }
    return 0;
}

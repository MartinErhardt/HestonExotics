/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include "WebAPI.h"
#include "../simdjson.h"
#include"BSM.h"
//#include <math.h>
using namespace WebInterface;

unsigned int days_to_date(const std::string& expiration_date);
size_t write_data(void* buffer, size_t size, size_t nmemb, WebInterface::JSON_buffer* buf){
    char* new_loc=NULL;
    if (!(new_loc= (char*)malloc(buf->size_buf+nmemb*size+1+simdjson::SIMDJSON_PADDING))) return 0;
    if (buf->buf) memcpy(new_loc,buf->buf,buf->size_buf);
    memcpy(new_loc+buf->size_buf,buffer,nmemb*size);
    *(new_loc+buf->size_buf+nmemb*size)='\0';
    if (buf->buf) free(buf->buf);
    buf->buf=new_loc;
    buf->size_buf=buf->size_buf+nmemb*size;
    return size * nmemb;
} //TODO number of workdays excluding holidays
unsigned int days_to_date(const std::string& expiration_date){
    struct std::tm then={};
    std::istringstream ss(expiration_date.c_str());
    ss>>std::get_time(&then, "%Y-%m-%d");
    unsigned int raw_diff=static_cast<unsigned int>(std::difftime(std::mktime(&then),time(0))/(60*60*24));
    unsigned int work_day_diff=5*static_cast<unsigned int>(raw_diff/7)+raw_diff%7;
    //std::cout<<"raw_diff"<<raw_diff<<"\t work_day_diff"<<work_day_diff<<'\n';
    return work_day_diff;
}
WebAPI::WebAPI(const std::string& access_code):token(access_code){
    buf.size_buf=0;
    buf.buf=NULL;
    header_list=NULL;
    simdjson::ondemand::parser JSONParser;
    if(!(curl = curl_easy_init())) throw std::runtime_error("could not initialize curl");
    header_list = curl_slist_append(header_list, AUTH_HEADER(token));
    header_list = curl_slist_append(header_list, JSON_HEADER);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header_list);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &buf);
}

void WebAPI::download_to_buf(const std::string& url){
    CURLcode res;
    free(buf.buf);
    buf.size_buf=0;
    buf.buf=NULL;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    if((res = curl_easy_perform(curl)) != CURLE_OK) throw std::runtime_error("curl failed in download_to_buf");
    //std::cout<<"WTF: "<<buf.buf<<'\n';
}
std::unique_ptr<std::list<std::string>> WebAPI::parse_expiries(){
    //std::cout<<buf.buf<< '\n';
    auto expiries = std::unique_ptr<std::list<std::string>>(new std::list<std::string>);
    auto doc = JSONParser.iterate(buf.buf, strlen(buf.buf), buf.size_buf+1+simdjson::SIMDJSON_PADDING);
    auto JSONList=doc.get_object()["expirations"];
    if(JSONList.type()==simdjson::ondemand::json_type::null) throw APIError(); // This happens if there have no valid underlying
    for(auto date : JSONList["date"]) expiries->push_back(std::string(std::string_view(date)));
    return expiries;
}
void WebAPI::parse_option_chain(options_chain& opt_chain,const char* vol_type, int vol_n){
    auto doc = JSONParser.iterate(buf.buf, strlen(buf.buf), buf.size_buf+1+simdjson::SIMDJSON_PADDING);
    for(auto opt : doc.get_object()["options"]["option"]){
        std::string option_type=std::string(std::string_view(opt["option_type"]));
        auto strike_obj=opt["strike"];
        auto price_obj=opt["ask"];
        int64_t current_vol=opt[vol_type].get_int64();
        if(option_type=="call"
            &&((opt["volume"]).type()!=simdjson::ondemand::json_type::null)
            &&(strike_obj.type()!=simdjson::ondemand::json_type::null)
            &&(price_obj.type()!=simdjson::ondemand::json_type::null)
            &&(current_vol>=vol_n)
        )
        {
            option * new_opt =new option();
            new_opt->volume=current_vol;
            new_opt->strike=static_cast<ffloat>(strike_obj.get_double());
            new_opt->price=static_cast<ffloat>(price_obj.get_double());
            new_opt->bid=static_cast<ffloat>(opt["bid"].get_double());
            opt_chain.options.push_back(*new_opt);
            opt_chain.min_strike=std::min(new_opt->strike, opt_chain.min_strike);
            opt_chain.max_strike=std::max(new_opt->strike, opt_chain.max_strike);
        }
    }
}
ffloat WebAPI::parse_stock_quote(){
    auto doc = JSONParser.iterate(buf.buf, strlen(buf.buf), buf.size_buf+1+simdjson::SIMDJSON_PADDING);
    auto obj= doc.get_object()["quotes"]["quote"]; //TODO Exception handling
    auto price_obj=obj["bid"];
    if(price_obj.type()!=simdjson::ondemand::json_type::null) return static_cast<ffloat>(price_obj.get_double());
    else throw APIError();
}
std::list<options_chain>* WebAPI::get_all_option_chains(const std::string& underlying,const char* vol_type, int vol_n){
    auto all_chains=new std::list<options_chain>(); 
    std::cout<<"Download expiries...";
    download_to_buf(EXP_URL(underlying));
    std::cout<<"Done\n";
    std::cout<<"Parse expiries...";
    auto expiries =parse_expiries();
    std::cout<<"Done\n";
    for(const std::string& date : *expiries){
        unsigned int cur_days_to_date=days_to_date(date);
        ffloat cur_time_to_date = static_cast<ffloat>(days_to_date(date))/trading_days;
        if (cur_time_to_date<=EXP_LB) continue;
        options_chain& new_opt_chain=*(new options_chain(cur_days_to_date,cur_time_to_date));
        std::cout<<"Download all options expiring on "<<date<<"...";
        download_to_buf(OPTIONS_CHAIN_URL(underlying,date));
        std::cout<<"Done\n";
        std::cout<<"Parse all options expiring on "<<date<<"...";
        parse_option_chain(new_opt_chain,vol_type,vol_n);
        std::cout<<"Done\n";
        all_chains->push_back(std::move(new_opt_chain));
    }
    return all_chains;
}
ffloat WebAPI::get_stock_quote(const std::string& stock){
    download_to_buf(QUOTE_URL(stock));
    return parse_stock_quote();
}
WebAPI::~WebAPI()
{
    curl_easy_cleanup(curl);
    curl_slist_free_all(header_list);
    free(buf.buf);
}

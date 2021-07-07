/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include "WebAPI.h"
#include "../simdjson.h"
//#include <math.h>
using namespace WebInterface;
//using namespace simdjson;
#define YMD_DATE_LENGTH 10
size_t write_data(void *buffer, size_t size, size_t nmemb, WebInterface::JSON_buffer*userp);
unsigned int days_to_date(const std::string& expiration_date);
size_t write_data(void *buffer, size_t size, size_t nmemb, WebInterface::JSON_buffer* buf)
{
    char* new_loc=NULL;
    if (!(new_loc= (char*)malloc(buf->size_buf+nmemb*size+1+simdjson::SIMDJSON_PADDING))) return 0;
    if (buf->buf) memcpy(new_loc,buf->buf,buf->size_buf);
    memcpy(new_loc+buf->size_buf,buffer,nmemb*size);
    *(new_loc+buf->size_buf+nmemb*size)='\0';
    if (buf->buf) free(buf->buf);
    buf->buf=new_loc;
    buf->size_buf=buf->size_buf+nmemb*size;
    //std::cout<<buf->buf<< '\n';
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
WebAPI::WebAPI(const std::string& access_code):token(access_code)
{
    std::cout << token << '\n';
    buf.size_buf=0;
    buf.buf=NULL;
    header_list=NULL;
    simdjson::ondemand::parser JSONParser;
    if(!(curl = curl_easy_init())) throw APIError();
    header_list = curl_slist_append(header_list, AUTH_HEADER(token));
    header_list = curl_slist_append(header_list, JSON_HEADER);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header_list);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &buf);
}

void WebAPI::download_to_buf(const std::string& url)
{
    CURLcode res;
    free(buf.buf);
    buf.size_buf=0;
    buf.buf=NULL;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    if((res = curl_easy_perform(curl)) != CURLE_OK) throw APIError();
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
void WebAPI::parse_option_chain(std::list<option>& options, unsigned int days_to_expire){
    //std::cout<<buf.buf<< '\n';
    auto doc = JSONParser.iterate(buf.buf, strlen(buf.buf), buf.size_buf+1+simdjson::SIMDJSON_PADDING);
    for(auto opt : doc.get_object()["options"]["option"]){
        int64_t vol;
        std::string option_type=std::string(std::string_view(opt["option_type"]));
        if(option_type=="call"
            &&!(opt["bid"].type()==simdjson::ondemand::json_type::null)
            &&!(opt["ask"].type() ==simdjson::ondemand::json_type::null)
            &&!(opt["strike"].type() ==simdjson::ondemand::json_type::null)
            &&!(opt["volume"].type() ==simdjson::ondemand::json_type::null)
            &&!(opt["last"].type() ==simdjson::ondemand::json_type::null)
            &&!((vol=opt["volume"].get_int64())==0))
        {
            option * new_opt=new option();
            new_opt->strike=static_cast<ffloat>(opt["strike"].get_double());
            new_opt->days_to_expiry=days_to_expire;
            ffloat bid=static_cast<ffloat>(opt["bid"].get_double());
            ffloat ask=static_cast<ffloat>(opt["ask"].get_double());
            new_opt->volume=vol;
            new_opt->price=static_cast<ffloat>(opt["last"].get_double());
            new_opt->ask=ask;
            new_opt->bid=bid;
            options.push_back(*new_opt);
            //std::cout<<"stock price: "<<S_1<<"\tprice: "<<new_opt->price<<"\tstrike: "<<new_opt->strike<<"\trisk free rate: "<<risk_free<<"\ttrading volume: "<<new_opt->volume<<"\tdays to expiry: "<<new_opt->days_to_expiry<<"\t spread: "<<new_opt->ask-new_opt->bid<<'\n';
        }
    }
}
ffloat WebAPI::parse_stock_quote(){
    auto doc = JSONParser.iterate(buf.buf, strlen(buf.buf), buf.size_buf+1+simdjson::SIMDJSON_PADDING);
    auto obj= doc.get_object()["quotes"]["quote"]["last"]; //TODO Exception handling
    if(!(obj.type()==simdjson::ondemand::json_type::null)) return static_cast<ffloat>(obj.get_double());
    else throw APIError();
    /*if(!(obj["bid"].type()==simdjson::ondemand::json_type::null)
            &&!(obj["ask"].type() ==simdjson::ondemand::json_type::null)){
        ffloat bid=static_cast<ffloat>(obj["bid"].get_double());
        ffloat ask=static_cast<ffloat>(obj["ask"].get_double());
        return (bid+ask)/2;
    } else throw APIError();*/
}
std::unique_ptr<std::list<option>> WebAPI::get_all_option_chains(const std::string& underlying){
    auto options=std::make_unique<std::list<option>>(); 
    std::cout<<"Download expiries...";
    download_to_buf(EXP_URL(underlying));
    std::cout<<"Done\n";
    std::cout<<"Parse expiries...";
    auto expiries =parse_expiries();
    std::cout<<"Done\n";
    unsigned int days_to_expire;
    for(const std::string& date : *expiries){
        std::cout<<"Download all options expiring on "<<date<<"...";
        download_to_buf(OPTIONS_CHAIN_URL(underlying,date));
        std::cout<<"Done\n";
        std::cout<<"Parse all options expiring on "<<date<<"...";
        days_to_expire=days_to_date(date);
        parse_option_chain(*options,days_to_expire);
        std::cout<<"Done\n";
    }
    return options;
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

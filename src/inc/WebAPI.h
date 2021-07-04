/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#ifndef WEB_API_H
#define WEB_API_H
#ifndef WIN32
#include <unistd.h>
#endif
#include <curl/curl.h>
#include "../../simdjson.h"
#include <list>
#include "Types.h"
namespace WebInterface
{
    #define EXP_URL(SYMB) std::string("https://sandbox.tradier.com/v1/markets/options/expirations?symbol=")+SYMB+std::string("&includeAllRoots=false&strikes=false")
    #define QUOTE_URL(SYMB) std::string("https://sandbox.tradier.com/v1/markets/quotes?symbols=")+SYMB
    #define OPTIONS_CHAIN_URL(SYMB,DATE) (std::string("https://sandbox.tradier.com/v1/markets/options/chains?symbol=") + SYMB + std::string("&expiration=") +DATE+ std::string("&greeks=false")).c_str()
    #define AUTH_HEADER(TOKEN) (std::string("Authorization: Bearer ") + TOKEN).c_str()
    #define JSON_HEADER "Accept: application/json"
    typedef struct{
        unsigned int days_to_expiry;
        ffloat price;
        ffloat strike;
        int64_t volume;
    } option;
    typedef struct{
        size_t size_buf;
        char * buf;
    } JSON_buffer;
    struct APIError : public std::exception
    {
        const char * what () const throw ()
        {
            return "WebAPI failed";
        }
    };
    
    class WebAPI{
        struct curl_slist *header_list;
        CURL *curl;
        JSON_buffer buf;
        const std::string& token;
        simdjson::ondemand::parser JSONParser;
        void download_to_buf(const std::string& url);
        std::unique_ptr<std::list<std::string>> parse_expiries();
        void parse_option_chain(std::list<option>& options, unsigned int days_to_date);
        ffloat parse_stock_quote();
        public:
            WebAPI(const std::string& access_code);
            std::unique_ptr<std::list<option>> get_all_option_chains(const std::string& underlying);
            ffloat get_stock_quote(const std::string& stock);
            ~WebAPI();
    };
}

#endif

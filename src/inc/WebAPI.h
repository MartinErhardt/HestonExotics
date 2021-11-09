/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once 
#ifndef WIN32
#include <unistd.h>
#endif
#include <curl/curl.h>
#include "../../simdjson.h"
#include <list>
#include "Types.h"

namespace WebInterface
{
    /** /addtogroup WebInterface
     * @{
     */
    #define EXP_URL(SYMB) std::string("https://sandbox.tradier.com/v1/markets/options/expirations?symbol=")+SYMB+std::string("&includeAllRoots=false&strikes=false") ///< URL for fething all expiry dates in tradier sandbox API 
    #define QUOTE_URL(SYMB) std::string("https://sandbox.tradier.com/v1/markets/quotes?symbols=")+SYMB ///< URL for stock quotes in tradier sandbox API
    #define OPTIONS_CHAIN_URL(SYMB,DATE) (std::string("https://sandbox.tradier.com/v1/markets/options/chains?symbol=") + SYMB + std::string("&expiration=") +DATE+ std::string("&greeks=false")).c_str() ///< URL for stock quotes in tradier sandbox API U
    #define AUTH_HEADER(TOKEN) (std::string("Authorization: Bearer ") + TOKEN).c_str() ///< authentication request in tradier sandbox API
    #define JSON_HEADER "Accept: application/json" ///< HTTP Header for for authentication in tradier sandbox API
    /**
     * @brief wrapper for the download buffer
     */
    typedef struct{
        size_t size_buf;    //<size of download buffer
        char * buf;         //<download buffer containing the tradier API responses in JSON format
    } JSON_buffer;
    /**
     * @brief very minimalistic Error class for handling tradier API errors.
     */
    struct APIError : public std::exception
    {
        const char * what () const throw ()
        {
            return "WebAPI failed";
        }
    };
    
    /**
    * @brief This is the WebAPI class, which is used to encapsulate all required resources and temporary data for downloading and parsing live market through libcurl, tradier WebAPI and simdjson.
    */
    class WebAPI{
        struct curl_slist *header_list;
        CURL *curl;
        JSON_buffer buf;
        const std::string& token;
        simdjson::ondemand::parser JSONParser;
        /** 
         * private method to obtain responses from tradier in JSON format through libcurl
         * @param url the API request URL
         * @sideeffect insert pointer to the API response in this.buf.buf and write size of response in this.buf.size_buf 
         * @throws std::runtime_error if download fails
         */
        void download_to_buf(const std::string& url);
        /**
         * private method to parse the list of expiries in the tradier JSON response to EXP_URL
         * 
         * only call, when tradier API response to EXP_URL can be found in this.download_to_buf
         * 
         * @return smart pointer to the list of all expiries given in the tradier API response in this.buf.buf with expiry dates in Y-m-d string format
         * @throws APIError if tradier API response does not satisfy JSON format corresponding to EXP_URL request
         */
        std::list<std::string> parse_expiries();
        /**
         * private method to parse the list of options in the tradier JSON response to OPTIONS_CHAIN_URL
         * 
         * only call, when tradier API response to OPTIONS_CHAIN_URL can be found in this.buf
         * 
         * @return smart pointer to the list of all expiries in string in Y-m-d format
         * @param opt_chain options chain in which to insert the information in buf.buf
         * @param vol_type either volume for trading volume or open_interest for trades outstanding
         * @param vol_n number of contracts specified vol_type
         * @sideeffect update opt_chain
         * @throws APIError if tradier API response does not satisfy JSON format corresponding to OPTIONS_CHAIN_URL request
         */
        void parse_option_chain(options_chain& opt_chain,const char* vol_type, int vol_n);
        /**
         * private method to parse the stock quote obtained through tradier API
         * 
         * only call, when tradier API response to QUOTE_URL can be found in this.buf
         * 
         * @return value of the stock in this.buf
         * @throws APIError if tradier API response does not satisfy JSON format corresponding to QUOTE_URL request
         */
        ffloat parse_stock_quote();
        public:
            WebAPI(const std::string& access_code);
            /**
             * public method to obtain all call options on a given underlying
             * 
             * Calls download_to_buf, then parse_option_chain on every single expiry_date  
             * @param underlying name of underlying
             * @param vol_type either volume for trading volume or open_interest for trades outstanding
             * @param vol_n number of contracts specified vol_type
             * @return list of options_chain returned by parse_option_chain
             * @throws std::runtime_error if curl can not be initialized
             */
            std::list<options_chain> get_all_option_chains(const std::string& underlying,const char* vol_type, int vol_n);
            /**
             * public method to obtain the value of any given stock listing. 
             * 
             * Calls download_to_buf, then parse_stock_quote 
             * @param stock name of the stock
             * @return value of parse_stock_quote
             */
            ffloat get_stock_quote(const std::string& stock);
            ~WebAPI();
    };
    /** @} 
     */
}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include"HDistribution.h"
#include<memory>
#include<list>
#include"SWIFT.h"
#include"BSM.h"


/**
 * @brief The objective of this struct type is to create a new data structure, that unlike options_chain also contains more implementation specific attributes like the distribution at expiry and the SWIFT method used for pricing
 */
struct ED{
    const options_chain& opts;  //<European Call options to price using the SWIFT method
    HDistribution* distr;       //<price distribution at the point of expiry of all the opts
    SWIFT* pricing_method;      //<SWIFT method used to price opts
    ED(const options_chain& opts,HDistribution* init_distr, SWIFT* init_pricing_method): opts(opts),distr(init_distr),pricing_method(init_pricing_method){}
    ~ED(){delete distr; delete pricing_method;}
};
typedef ED expiry_data;

/**
 * @brief adata struct type used to store additional data passed through the dlevmar calls
 */ 
struct AS{
    ffloat S;                           //<price of underlying
    ffloat* real_prices;                //<pointer to C-array of prices observed in the market(initialization optional)
    std::list<expiry_data>& exp_list;   //<list of expiry data. Use of reference to allow overloaded [] opertor with resource owned by AS 
    ~AS(){delete &exp_list;}            //moved when push_back 
};
typedef AS adata_s;
/**
 * custom shift operator used for debugging in the standard output, output format:
 * S: (stock price) [tab]emp price: (price observed in the market)[tab]ask: (ask price of call option)[tab]bid: (bid price of option)[tab]strike: (call strike)[tab]volume: (volume open interest or trading volume on last/current day depending on program aguments)[tab]days left: (days left) 
 */
std::ostream& operator<<(std::ostream& out, adata_s const& as);
/**
 * @brief function determining the prices of European Call options in the Heston model using the SWIFT method
 * 
 * Right now this function (or the rather local function update_adata) deletes and reinitializes the pricing_method in expiry_data each time, the model parameters change. This is done to achieve higher precision, because the integral bounds lower and upper depend on the model parameters. That being said it is perfectly fine to reuse the coefficients U_mj by not deleting and reinitializing a new SWIFT method every time the parameters change(and also better in terms of memory fragmentation).
 * Note that this C-function is passed to dlevmar. 
 * 
 * @param p pointer to the current parameters
 * @param x pointer to the output price data
 * @param m number of parameters to calibrate (in this case =1) 
 * @param n_observations number of options to price
 * @param adata pointer to the call invariant data
 * @throws std::runtime_error if SWIFT::price_opts does not price all options
 * @sideeffect write results x
 */
void get_prices_for_levmar(ffloat *p, ffloat *x, int m, int n_observations, void * adata);
/**
 * @brief function determining the price gradient of European Call options in the Heston model using the SWIFT method
 * 
 * Right now this function (or the rather local function update_adata) deletes and reinitializes the pricing_method in expiry_data each time, the model parameters change. This is done to achieve higher precision, because the integral bounds lower and upper depend on the model parameters. That being said it is perfectly fine to reuse the coefficients U_mj by not deleting and reinitializing a new SWIFT method every time the parameters change(and also better in terms of memory fragmentation).
 * Note that this C-function passed to dlevmar.
 * 
 * @param p pointer to the current parameters
 * @param jac pointer to the output price gradient data
 * @param m number of parameters to calibrate (in this case =5) 
 * @param n_observations number of options to obtain price gradient from
 * @param adata pointer to the call invariant data
 * @throws std::runtime_error if SWIFT::price_opts_grad does not compute the gradient information for all options
 * @sideeffect write results to jac
 */
void get_jacobian_for_levmar(ffloat *p, ffloat *jac, int m, int n_observations, void * adata);
/**
 * @brief function calibrating the Heston model to a stock given its price and various European call options to different expiry dates on that stock
 * 
 * This function first initializes an adata_s data structure based on market_data. Then it calls the Levenberg-Marquardt routine with that adata, get_prices_for_levmar and get_jacobian_for_levmar. The algorithm stops if the gradient is lesser equal than 1e^-10. 
 * 
 * Lastly note that by the original paper from Yiran Li, it is noted, that pricing european calls in the Heston model is numerically ill conditioned in the parameter v_m, when compared to the parameter sigma, but in particular compared to kappa, whose condition number is diminished by a factor of 1e-6. Conversely for the calibration model, this means that calibration is numerically ill-conditioned in the parameter kappa frequently resulting in large values of kappa>15  .
 * 
 * @param S stock price
 * @param market_data various options chains containing European call options to different expiry dates
 * @return smart pointer to new heston model parameters
 * @throws std::runtime_error if levmar fails for various reasons
 */
HParams calibrate(const ffloat S,const std::list<options_chain>& market_data);

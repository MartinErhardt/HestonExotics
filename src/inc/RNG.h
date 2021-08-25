/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<stdint.h>
#include"Types.h"
#define SEEDTYPE uint64_t
#define ALIGN (1<<7) ///<heap alignment
#include"shishua.h"
/**
 * @brief simple wrapper class for the shishua pseudo rng
 * 
 * Internally the wrapper is implemented through a very simple ring puffer
 *
 * @see  https://github.com/espadrine/shishua
 */
class RNG{
    ffloat * buf_start_u;
    ffloat * buf_end_u;
    ffloat * buf_cur_u;
    ffloat * buf_start_g;
    ffloat * buf_end_g;
    ffloat * buf_cur_g;
    SEEDTYPE my_seed[4]={0};
    prng_state s;
    ffloat* setup_u(ffloat* buf_start_setup,ffloat* buf_end_setup);
    ffloat* setup_g();
    public:
        /**
         * Constructor 
         * @param size size of ring puffer
         * @param seed seed for shishua
         */
        RNG(size_t size,unsigned int seed);
        /**
         * @brief Computes a standard normally distributed random variable
         * @return standard normally distributed random variable
         */
        ffloat get_grand(){
            if(buf_cur_g==buf_end_g) buf_cur_g=setup_g();
            return *(buf_cur_g++);
        }
        /**
         * @brief Computes a random variable that is distributed uniformly on the unit interval
         * @return random variable that is distributed uniformly on the unit interval
         */
        ffloat get_urand(){
            if(buf_cur_u==buf_end_u) buf_cur_u=setup_u(buf_start_u,buf_end_u);
            return *(buf_cur_u++);
        }
        ~RNG();
};

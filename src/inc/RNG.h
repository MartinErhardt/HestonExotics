/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<stdint.h>
#include"Types.h"
#define SEEDTYPE unsigned int
#define ALIGN (1<<7)
#include"shishua.h"

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
        RNG(size_t size,unsigned int seed);
        ffloat get_grand();
        ffloat get_urand();
        ~RNG();
};

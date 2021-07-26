/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<stdint.h>
#include"Types.h"
#define SEEDTYPE unsigned int
#define ALIGN (1<<6)
#include"shishua.h"

class RNG{
    ffloat * buf_start;
    ffloat * buf_end;
    ffloat * buf_cur;
    SEEDTYPE my_seed[4]={0};
    prng_state s;
    void setup();
    public:
        RNG(size_t size,unsigned int seed);
        ffloat get_grand();
        ffloat get_urand();
        ~RNG();
};

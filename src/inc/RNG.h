/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<stdint.h>
#include"Types.h"
class RNG{
    ffloat * buf_start;
    ffloat * buf_end;
    public:
        ffloat * buf_cur;
        RNG(size_t size);
        void update();
};

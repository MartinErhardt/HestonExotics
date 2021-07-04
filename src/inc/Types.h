/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#ifndef TYPES_H
#define TYPES_H
#include<stdint.h>
typedef double ffloat;
typedef struct{
    unsigned int days_to_expiry;
    ffloat price;
    ffloat strike;
    int64_t volume;
} option;
#endif

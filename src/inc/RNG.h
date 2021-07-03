/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#ifndef RNG_H
#define RNG_H

class RNG{
    double * buf_start;
    double * buf_end;
    public:
        double * buf_cur;
        RNG(size_t size);
        void update();
};
#endif

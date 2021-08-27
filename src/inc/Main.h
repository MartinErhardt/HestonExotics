/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
enum program_option{
    CALIBRATE=0,
    PRICE=1,
    DOWNLOAD=2,
    TEST=3,
    HELP=4
};
const char help_info[]="hexo [{-c <underlying> [<volume type> <# volume>]}...] [{-d <underlying> [<volume type> <# volume>]}...] [-h] [{-p <underlying> <option type> {all|<days to expiration> <strike>} }...] [{-t <test>}...]";
struct SynapsisError : public std::exception
{
    const char * what () const throw ()
    {
        return help_info;
    }
};

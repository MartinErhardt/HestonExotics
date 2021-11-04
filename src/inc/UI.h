/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<string>
namespace UI{
    #define PROG_BAR_WIDTH 100
    #define PROG_CHAR '#'

    class ProgressBar{
        char prog_bar[PROG_BAR_WIDTH+1];
        std::string message;
        unsigned int final_value;
        unsigned int current_value;
    public:
        ProgressBar(std::string&& message_i,unsigned int final_value_i);
        void update(unsigned int i);
        void finish();
    };
}

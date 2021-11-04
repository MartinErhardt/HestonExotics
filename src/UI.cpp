/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include<cstdio>
#include"UI.h"
using namespace UI;
ProgressBar::ProgressBar(std::string&& message_i,unsigned int final_value_i){
    message=std::move(message_i);
    for(unsigned int i=0;i<PROG_BAR_WIDTH;i++) prog_bar[i]=PROG_CHAR;
    prog_bar[PROG_BAR_WIDTH]='\0';
    current_value=0;
    final_value=final_value_i;
}
void ProgressBar::update(unsigned int i){
    double prog= ((double)(current_value=current_value+i))/((double)final_value);
    int visible_prog=(int)(prog*PROG_BAR_WIDTH);
    printf("\r%s\t[%.*s%*s] %3d%%\t(%d,%d)", message.c_str(),visible_prog, prog_bar, PROG_BAR_WIDTH-visible_prog, " ",(int) (prog*100),current_value,final_value);
    fflush(stdout);
}
void ProgressBar::finish(){
    printf("\r%s\t[%s] 100%%\t(%d,%d)\n",message.c_str(),prog_bar,current_value,final_value);
}

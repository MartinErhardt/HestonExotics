/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"RNG.h"
//#include"../as241.f"
//extern "C" prng_state prng_init(SEEDTYPE seed[4]);
extern "C" double ppnd16(double *,int*);

RNG::RNG(size_t size,unsigned int seed){
    my_seed[0]=seed;
    if((size*sizeof(ffloat))&(ALIGN-1)) throw std::runtime_error("RNG: init size does not satisfy alignment"); 
    buf_start=(ffloat*)aligned_alloc(ALIGN,size*sizeof(ffloat));
    buf_end=buf_start+size;
    s = prng_init(my_seed);
    setup();
}
void RNG::setup(){
    buf_cur=buf_start;
    prng_gen(&s, (uint8_t*) buf_start, sizeof(ffloat)*(buf_end-buf_start));
    ffloat* x2=buf_start;
    do *x2=((ffloat)((fuint)*x2))/((ffloat)((fuint)-1)); while(++x2<buf_end); //TODO too many casts
}
ffloat RNG::get_urand(){
    if(buf_cur==buf_end) setup();
    return *(buf_cur++);
}
ffloat RNG::get_grand(){
    int ifault=0;
    if(buf_cur==buf_end) setup();
    double retval=ppnd16(buf_cur++,&ifault);
    if(ifault) throw std::runtime_error("Quantile function failed!");
    return retval;
}
RNG::~RNG(){
    free(buf_start);
}

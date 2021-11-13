/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"RNG.h"

extern "C" double ppnd16(double *,int*);

RNG::RNG(size_t size,unsigned int seed){
    my_seed[0]=seed;
    assert(!((size*sizeof(ffloat))&(ALIGN-1)));
/**
#ifdef __MINGW64__
    buf_start_g=(ffloat*)_aligned_malloc(size*sizeof(ffloat),ALIGN);    
	buf_start_u=(ffloat*)_aligned_malloc(size*sizeof(ffloat),ALIGN);
#elif defined(__GNUC__)
    buf_start_g=(ffloat*)aligned_alloc(ALIGN,size*sizeof(ffloat));    
	buf_start_u=(ffloat*)aligned_alloc(ALIGN,size*sizeof(ffloat));
#endif
    **/
    buf_start_g=new(std::align_val_t{ALIGN}) ffloat[size];
    buf_start_u=new(std::align_val_t{ALIGN}) ffloat[size];
    buf_end_u=buf_start_u+size;
    buf_end_g=buf_start_g+size;
    prng_init(&s, my_seed);
    buf_cur_u=setup_u(buf_start_u,buf_end_u);
    buf_cur_g=setup_g();
}
ffloat* RNG::setup_u(ffloat* buf_start_setup,ffloat* buf_end_setup){
    prng_gen(&s, (uint8_t*) buf_start_setup, sizeof(ffloat)*(buf_end_setup-buf_start_setup));
    ffloat* x2=buf_start_setup;
    do *x2=((ffloat)*((fuint*)x2))/((ffloat)((fuint)-1)); while(++x2<buf_end_setup);
    return buf_start_setup;
}
ffloat* RNG::setup_g(){
    setup_u(buf_start_g,buf_end_g);
    buf_cur_g=buf_start_g;
    ffloat* x2=buf_start_g;
    do{ int ifault=0;
        *x2=ppnd16(x2,&ifault);
        assert(!ifault);
    }while((++x2)!=buf_end_g);
    return buf_start_g;
}
RNG::~RNG(){
    operator delete[](buf_start_g, std::align_val_t{ALIGN});
    operator delete[](buf_start_u, std::align_val_t{ALIGN});
}

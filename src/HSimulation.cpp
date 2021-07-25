/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"HDistribution.h"
#include"HSimulation.h"
#define PSI_C 1.5
template <typename SchemeParams,class Scheme>
SDE<2>& HQEAnderson<SchemeParams,Scheme>::operator++() {
    auto const& [v_0,theta,rho,kappa,eps] = params;
    ffloat delta=step_width(&state);
    ffloat gamma_1=.5;
    ffloat gamma_2=.5;
    ffloat discount=std::exp(-kappa*delta);
    ffloat m=theta+(state.cur[1]-theta)*discount;
    ffloat sp2=state.cur[1]*eps*eps*discount/kappa*(1-discount)+theta*eps*eps/(2*kappa)*discount*discount;
    ffloat Psi=sp2/(m*m);
    ffloat p=(Psi-1)/(Psi+1);
    ffloat beta=2/(m*(Psi+1));
    ffloat bp2=2/Psi-1+std::sqrt(2/Psi*(2/Psi-1));
    ffloat b=std::sqrt(bp2);
    ffloat a=m/(1+bp2);
    //auto Psi_inv = [](ffloat u, ffloat p, ffloat beta){
    //        return 
    //};
    state.prev[1]=state.cur[1];
    if(Psi>PSI_C){
        ffloat Z_V=rng->get_grand();
        state.cur[1]=a*(b+Z_V)*(b+Z_V);
    } else{
        ffloat U_V=rng->get_urand();
        state.cur[1]=p<U_V? std::log((1-p)/(1-U_V))/beta:0.;
    }
    ffloat K_0=-rho*kappa*theta/eps*delta;
    ffloat K_1=gamma_1*delta*(kappa*rho/eps-.5)-rho/eps;
    ffloat K_2=gamma_2*delta*(kappa*rho/eps-.5)+rho/eps;
    ffloat K_3=gamma_1*delta*(1-rho*rho);
    ffloat K_4=gamma_2*delta*(1-rho*rho);
    log_X=log_X+K_0+K_1*state.prev[1]+K_2*state.cur[1]+std::sqrt(K_3*state.prev[1]+K_4*state.cur[1])*rng->get_grand();
    state.prev[0]=state.cur[0];
    state.cur[0]=std::exp(log_X);
    state.time+=delta;
    return *this;
}

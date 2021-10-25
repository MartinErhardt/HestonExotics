/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#if __INCLUDE_LEVEL__
  #error "Don't include this file"
#endif

#include"HSimulation.cpp"
#include"HSimulation.h"
#include"AsianContract.h"

template class HSimulation::HQEAnderson<ffloat,AAsianCallNonAdaptive>;
template class HSimulation::PricingTool<HSimulation::HQEAnderson<ffloat,AAsianCallNonAdaptive>>;

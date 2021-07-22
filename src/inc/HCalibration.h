#pragma once 
#include"HDistribution.h"
#include<memory>
#include<list>

void pricing_test();
void gradient_test();
void levmar_test();
std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data);

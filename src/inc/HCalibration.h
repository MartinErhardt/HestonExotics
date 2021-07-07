#pragma once 
#include"HDistribution.h"
#include<memory>
#include<list>

std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<options_chain>& market_data);

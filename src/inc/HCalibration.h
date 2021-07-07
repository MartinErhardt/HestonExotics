#pragma once 
#include"HDistribution.h"
#include<memory>
#include<list>

std::unique_ptr<HParams> calibrate(const ffloat S,const std::list<option>& market_data);

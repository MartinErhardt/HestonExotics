/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<vector>
enum tests{
    DISTRIBUTION=0,     //<distribution test
    PRICING=1,          //<pricing test
    GRADIENT=2,         //<gradient test
    LEVMAR=3,           //<levmar test
    RNG=4               //<rng test
};
/**arguments passed to -t option in order of tests*/
const std::vector<std::string> arguments={"distr","pricing","levmar","RNG"};
/**
 * @brief tests the computation of the characteristic function and its partial derivatives
 * @sideeffect stdout
 * @throws abort if test fails
 */
void distr_test();
/**
 * @brief tests the price computation using the SWIFT method
 * @sideeffect stdout
 * @throws abort if test fails
 */
void pricing_test();
/**
 * @brief tests the computation of the price gradients using the SWIFT method
 * @sideeffect stdout
 * @throws abort if test fails
 */
void gradient_test();
/**
 * @brief tests the Heston model calibration using levmar
 * @sideeffect stdout
 * @throws abort if test fails
 */
void levmar_test();
/**
 * @brief very simple test used to ensure correctness of the shishua wrapper.
 * It computes the expectation and sample variances of both uniform distributions on the unit interval and standard normal distributions using a Monte Carlo method up to 1e-4 accuracy.
 * The number of trial is around 4 billion. 
 * 
 * A more advanced test could consist of computing the rate of convergence of the Monte Carlo method. Since shishua suffices much harder statistical tests, this is not necesary.
 * @sideeffect stdout
 * @throws abort if test fails
 */
void rng_test();

/**addresses of the test_functions in order of tests*/
void* test_functions[]={(void*)&distr_test,(void*)&pricing_test,(void*)&levmar_test,(void*)&rng_test};

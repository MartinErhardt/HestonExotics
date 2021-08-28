/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include<map>

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
/**
 * @brief tests the quadratic exponential (QE) path simulation
 * @sideeffect stdout
 * @throws abort if test fails
 */
void simulation_test();
typedef void (*test_func)();

const std::map<std::string, test_func> test_eval {
        {"distr", distr_test},
        {"pricing", pricing_test},
        {"gradient",gradient_test},
        {"levmar", levmar_test},
        {"rng", rng_test},
        {"simulation",simulation_test}
};

#ifndef __TEST_GENERATOR_H
#define __TEST_GENERATOR_H

#define _USE_MATH_DEFINES
#include <cmath>
#include "Descriptor.h"

void randomFlightState(
    FLIGHT_STATE* flightState, 
    double_t xUpperBound = 1000, 
    double_t yUpperBound = 1000, 
    double_t degUpperBound = 2 * M_PI,
    double_t xLowerBound = 0, 
    double_t yLowerBound = 0,
    double_t degLowerBound = 0);

void randomTask(
    TASK* task, 
    double_t xUpperBound = 1000, 
    double_t yUpperBound = 1000, 
    double_t degUpperBound = 2 * M_PI,
    double_t xLowerBound = 0, 
    double_t yLowerBound = 0,
    double_t degLowerBound = 0);

#endif
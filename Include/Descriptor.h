#ifndef __DESCRIPTOR__H
#define __DESCRIPTOR__H

#include <cstdint>
#include <cmath>

typedef struct tTASK_PARAMETER {
    double_t deg;
    double_t x;
    double_t y;
} TASK_PARAMETER;
typedef struct tTASK_QUEUE {
    uint32_t num;
    TASK_PARAMETER *tasks;
} TASK_QUEUE;
typedef struct tTASK {
    uint32_t queueNum;
    uint32_t totalNum;
    TASK_QUEUE *taskQueue;
} TASK;
typedef struct tSINGLE_FLIGHT_STATE {
    double_t deg;
    double_t x;
    double_t y;
} SINGLE_FLIGHT_STATE;
typedef struct tFLIGHT_STATE {
    uint32_t num;
    SINGLE_FLIGHT_STATE *flightState;
} FLIGHT_STATE;

#endif
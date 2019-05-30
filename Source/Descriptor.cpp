#include "Descriptor.h"

void createTask(TASK* task, uint32_t queueNum, const uint32_t queueLen) {
    task->queueNum = queueNum;
    task->totalNum = 0;
    task->taskQueue = new TASK_QUEUE[queueNum];
    for(uint32_t queueNo = 0; queueNo < queueNum; queueNo++) {
        task->totalNum += queueLen;
        task->taskQueue[queueNo].num = queueLen;
        task->taskQueue[queueNo].tasks = new TASK_PARAMETER[queueLen];
    }
}

void createTask(TASK* task, uint32_t queueNum, const uint32_t* queueLen) {
    task->queueNum = queueNum;
    task->totalNum = 0;
    task->taskQueue = new TASK_QUEUE[queueNum];
    for(uint32_t queueNo = 0; queueNo < queueNum; queueNo++) {
        task->totalNum += queueLen[queueNo];
        task->taskQueue[queueNo].num = queueLen[queueNo];
        task->taskQueue[queueNo].tasks = new TASK_PARAMETER[queueLen[queueNo]];
    }
}

void destructTask(TASK* task) {
    for(uint32_t queueNo = 0; queueNo < task->queueNum; queueNo++) {
        delete [] task->taskQueue[queueNo].tasks;
    }
    delete [] task->taskQueue;
    task->queueNum = 0;
    task->totalNum = 0;
}

void createFlight(FLIGHT_STATE* flightState, uint32_t flightNum) {
    flightState->num = flightNum;
    flightState->flightState = new SINGLE_FLIGHT_STATE[flightNum];
}

void destructFlight(FLIGHT_STATE* flightState) {
    flightState->num = 0;
    delete [] flightState->flightState;
}
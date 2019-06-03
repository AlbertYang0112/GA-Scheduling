#ifndef GAHW_DEBUG_H
#define GAHW_DEBUG_H

#include <cstdio>

#ifdef __DEBUG__
#define DEBUG(format,...) printf("File: " __FILE__ ", Line: %05d: " format "\n", __LINE__, ##__VA_ARGS__)
#define DEBUG_BRIEF(format,...) printf("" format "", ##__VA_ARGS__)
#else
#define DEBUG(format,...)
#define DEBUG_BRIEF(format,...)
#endif

#endif //GAHW_DEBUG_H

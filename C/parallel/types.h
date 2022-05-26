#ifndef __TYPES_H__
#define __TYPES_H__

#include <stdint.h>

typedef uint16_t dimensionIndex_t;
typedef uint64_t pointIndex_t;

typedef struct {
    double self;
    double recv;
    double best;
} MinDistance;

typedef struct {
    pointIndex_t self[2];
    pointIndex_t recv[2];
} BestPair;

#endif
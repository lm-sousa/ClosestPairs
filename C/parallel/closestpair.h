#ifndef __CLOSESTPAIR_H__
#define __CLOSESTPAIR_H__

#include <float.h>

#include "types.h"

void closestPair(double *points, pointIndex_t bestPair[2], double *minDistance,
                 pointIndex_t numberOfPoints, dimensionIndex_t dims);

#endif
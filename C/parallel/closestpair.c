#include "closestpair.h"

#define getPointIndex(pointNumber) (dims * (pointNumber))

void closestPair(double *points, pointIndex_t bestPair[2], double *minDistance,
                 pointIndex_t numberOfPoints, dimensionIndex_t dims) {

    (*minDistance) = DBL_MAX;

    for (pointIndex_t anchorPoint = 0; anchorPoint < numberOfPoints;
         anchorPoint++) {
        for (pointIndex_t point = anchorPoint + 1; point < numberOfPoints;
             point++) {
            double pointDistance = 0.0;
            for (dimensionIndex_t dim = 0; dim < dims; dim++) {
                double dimensionDistance =
                    (points[getPointIndex(anchorPoint) + dim] -
                     points[getPointIndex(point) + dim]);
                pointDistance += dimensionDistance * dimensionDistance;
            }

            if (pointDistance < (*minDistance)) {
                bestPair[0] = anchorPoint;
                bestPair[1] = point;
                (*minDistance) = pointDistance;
            }
        }
    }
}

void closestPairMerge(double *a) { return; }
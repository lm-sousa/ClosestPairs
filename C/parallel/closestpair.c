#include "closestpair.h"

#define getPointIndex(pointNumber) (dims * (pointNumber))

double distanceCalculation(double *points, pointIndex_t pointA,
                           pointIndex_t pointB, dimensionIndex_t dims) {
    double pointDistance = 0.0;
    for (dimensionIndex_t dim = 0; dim < dims; dim++) {
        double dimensionDistance = (points[getPointIndex(pointA) + dim] -
                                    points[getPointIndex(pointB) + dim]);
        pointDistance += dimensionDistance * dimensionDistance;
    }

    return pointDistance;
}

void closestPair2(double *points, pointIndex_t bestPair[2], double *minDistance,
                  pointIndex_t numberOfPoints, dimensionIndex_t dims) {

    (*minDistance) = DBL_MAX;

    for (pointIndex_t anchorPoint = 0; anchorPoint < numberOfPoints;
         anchorPoint++) {
        for (pointIndex_t point = anchorPoint + 1; point < numberOfPoints;
             point++) {
            double pointDistance =
                distanceCalculation(points, anchorPoint, point, dims);

            if (pointDistance < (*minDistance)) {
                bestPair[0] = anchorPoint;
                bestPair[1] = point;
                (*minDistance) = pointDistance;
            }
        }
    }
}

void closestPair(double *points, pointIndex_t bestPair[2], double *minDistance,
                 pointIndex_t numberOfPoints, dimensionIndex_t dims) {

    if (numberOfPoints <= 3) {
        closestPair2(points, bestPair, minDistance, numberOfPoints, dims);
        return;
    }

    pointIndex_t leftBestPair[2];
    double leftMinDistance;

    pointIndex_t rightBestPair[2];
    double rightMinDistance;

    pointIndex_t stripBestPair[2];
    double stripMinDistance;

    pointIndex_t middlePoint = numberOfPoints / 2;

    closestPair(&points[0], leftBestPair, &leftMinDistance, middlePoint, dims);

    closestPair(&points[getPointIndex(middlePoint)], rightBestPair,
                &rightMinDistance, numberOfPoints - middlePoint, dims);

    if (leftMinDistance <= rightMinDistance) {
        (*minDistance) = leftMinDistance;
        bestPair[0] = leftBestPair[0];
        bestPair[1] = leftBestPair[1];
    } else {
        (*minDistance) = rightMinDistance;
        bestPair[0] = rightBestPair[0] + middlePoint;
        bestPair[1] = rightBestPair[1] + middlePoint;
    }

    pointIndex_t leftPointIndex = getPointIndex(middlePoint - 1);
    pointIndex_t middlePointIndex = leftPointIndex;
    pointIndex_t rightPointIndex = getPointIndex(middlePoint);
    double realMinDistance = sqrt(*minDistance);

    while (leftPointIndex > 0 && ((points[middlePointIndex] -
                                   points[leftPointIndex]) < realMinDistance)) {
        leftPointIndex -= dims;
    }

    while (rightPointIndex < (getPointIndex(numberOfPoints - 1)) &&
           ((points[rightPointIndex] - points[middlePointIndex]) <
            realMinDistance)) {
        rightPointIndex += dims;
    }

    closestPair2(&points[leftPointIndex], stripBestPair, &stripMinDistance,
                 (rightPointIndex - leftPointIndex) / dims + 1, dims);

    if (stripMinDistance < (*minDistance)) {
        (*minDistance) = stripMinDistance;
        bestPair[0] = stripBestPair[0] + leftPointIndex / dims;
        bestPair[1] = stripBestPair[1] + leftPointIndex / dims;
    }
}

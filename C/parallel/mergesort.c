#include "mergesort.h"

#define getPointIndex(pointNumber) (dimensions * (pointNumber))

void mergeSortRecursive(double *outputArray, double *inputArray,
                        pointIndex_t leftPoint, pointIndex_t rightPoint,
                        dimensionIndex_t dimensions) {

    if (leftPoint < rightPoint) {
        pointIndex_t middlePoint = (leftPoint + rightPoint) / 2;
        mergeSortRecursive(outputArray, inputArray, leftPoint, middlePoint,
                           dimensions);
        mergeSortRecursive(outputArray, inputArray, middlePoint + 1, rightPoint,
                           dimensions);
        merge(outputArray, inputArray, leftPoint, middlePoint, rightPoint,
              dimensions);
    }

    return;
}

void mergeSort(double *outputArray, double *inputArray,
               pointIndex_t inputArraySize, dimensionIndex_t dimensions) {

    pointIndex_t numberOfPoints = inputArraySize / dimensions;
    mergeSortRecursive(outputArray, inputArray, 0, numberOfPoints - 1,
                       dimensions);
    return;
}

void merge(double *outputArray, double *inputArray, pointIndex_t leftPoint,
           pointIndex_t middlePoint, pointIndex_t rightPoint,
           dimensionIndex_t dimensions) {

    pointIndex_t leftCursor = getPointIndex(leftPoint);
    pointIndex_t rightCursor = getPointIndex(middlePoint + 1);
    pointIndex_t outputCursor = getPointIndex(leftPoint);

    while ((leftCursor <= getPointIndex(middlePoint)) &&
           (rightCursor <= getPointIndex(rightPoint))) {
        if (outputArray[leftCursor] < outputArray[rightCursor]) {
            for (dimensionIndex_t i = 0; i < dimensions; i++) {
                inputArray[outputCursor++] = outputArray[leftCursor++];
            }
        } else {
            for (dimensionIndex_t i = 0; i < dimensions; i++) {
                inputArray[outputCursor++] = outputArray[rightCursor++];
            }
        }
    }

    while (leftCursor <= getPointIndex(middlePoint)) {
        for (dimensionIndex_t i = 0; i < dimensions; i++) {
            inputArray[outputCursor++] = outputArray[leftCursor++];
        }
    }

    while (rightCursor <= getPointIndex(rightPoint)) {
        for (dimensionIndex_t i = 0; i < dimensions; i++) {
            inputArray[outputCursor++] = outputArray[rightCursor++];
        }
    }

    for (pointIndex_t j = getPointIndex(leftPoint);
         j <= getPointIndex(rightPoint); j += dimensions) {
        for (dimensionIndex_t i = 0; i < dimensions; i++) {
            outputArray[j + i] = inputArray[j + i];
        }
    }

    return;
}
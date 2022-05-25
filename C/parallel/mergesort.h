#ifndef __MERGESORT_H__
#define __MERGESORT_H__

#include "types.h"

void mergeSort(double *outputArray, double *inputArray,
               pointIndex_t inputArraySize, dimensionIndex_t dimensions,
               dimensionIndex_t dimensionToSortBy);

void merge(double *outputArray, double *inputArray, pointIndex_t leftPoint,
           pointIndex_t middlePoint, pointIndex_t rightPoint,
           dimensionIndex_t dimensions, dimensionIndex_t dimensionToSortBy);

#endif
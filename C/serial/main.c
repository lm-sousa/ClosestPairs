#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define dimensionIndex_t uint16_t
#define pointIndex_t uint64_t

dimensionIndex_t dims = 2;
pointIndex_t numberOfPoints = 0;

double *dimensionData = NULL;
double *points = NULL;

#define dimensionMin(dimension) (dimensionData[2 * dimension + 0])
#define dimensionMax(dimension) (dimensionData[2 * dimension + 1])

#define getPointIndex(pointNumber) (dims * pointNumber)

clock_t timePoint[4];
double timeDiff;

int main(int argc, char const *argv[]) {

    timePoint[0] = clock();

    if (argc < 3) {
        printf("Usage: <executable> inputFile.txt outputFile.txt\n");
        return 1;
    }

    int numValuesRead;

    FILE *fi = fopen(argv[1], "r");
    if (fi == NULL) {
        printf("Error opening input file.\n");
        return -1;
    }

    FILE *fo = fopen(argv[2], "w");
    if (fi == NULL) {
        printf("Error opening output file.\n");
        fclose(fi);
        return -1;
    }

    numValuesRead = fscanf(fi, "%lu", &numberOfPoints);
    if (numValuesRead != 1) {
        fclose(fi);
        fclose(fo);
        printf("Error reading number of points from file.\n");
        return -1;
    }

    if (numberOfPoints < 2) {
        fclose(fi);
        fclose(fo);
        printf("Number of points needs to be at least 2.\n");
        return -1;
    }

    numValuesRead = fscanf(fi, "%hu", &dims);
    if (numValuesRead != 1) {
        fclose(fi);
        fclose(fo);
        printf("Error reading dimensions from file.\n");
        return -1;
    }

    dimensionData = (double *)malloc(2 * sizeof(double) * dims);
    if (dimensionData == NULL) {
        fclose(fi);
        fclose(fo);
        printf("Error allocating memory for dimension data.\n");
        return -1;
    }

    double *points = (double *)malloc(sizeof(double) * dims * numberOfPoints);
    if (points == NULL) {
        fclose(fi);
        fclose(fo);
        free(dimensionData);
        printf("Error allocating memory for points.\n");
        return -1;
    }

    for (dimensionIndex_t i = 0; i < dims; i++) {
        numValuesRead =
            fscanf(fi, "%lf, %lf", &dimensionMin(i), &dimensionMax(i));
        if (numValuesRead != 2) {
            fclose(fi);
            fclose(fo);
            free(dimensionData);
            free(points);
            printf("Error reading dimension boundaries from file.\n");
            return -1;
        }
    }

    for (pointIndex_t point = 0; point < numberOfPoints; point++) {
        for (dimensionIndex_t dim = 0; dim < dims; dim++) {

            numValuesRead =
                fscanf(fi, "%lE, ", &points[getPointIndex(point) + dim]);
            if (numValuesRead != 1) {
                fclose(fi);
                fclose(fo);
                free(dimensionData);
                free(points);
                printf("Error reading points from file.\n");
                return -1;
            }
        }
    }

    timePoint[1] = clock();

    pointIndex_t bestPair[2];
    double minDistance = DBL_MAX;

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

            if (pointDistance < minDistance) {
                bestPair[0] = anchorPoint;
                bestPair[1] = point;
                minDistance = pointDistance;
            }
        }
    }

    timePoint[2] = clock();

    for (uint8_t i = 0; i < 2; i++) {
        for (dimensionIndex_t dim = 0; dim < dims; dim++) {
            fprintf(fo, "%lf ", points[getPointIndex(bestPair[i]) + dim]);
        }
        fprintf(fo, "\n");
    }

    fclose(fi);
    fclose(fo);
    free(dimensionData);
    free(points);

    timePoint[3] = clock();

    printf("\nD: %hu; N: %lu\n", dims, numberOfPoints);
    printf("\nTime elapsed in miliseconds:\n");
    printf("\tReading files ----- : %6.3lf\n",
           (double)(timePoint[1] - timePoint[0]) / CLOCKS_PER_SEC * 1000);
    printf("\tFinding closest pair: %6.3lf\n",
           (double)(timePoint[2] - timePoint[1]) / CLOCKS_PER_SEC * 1000);
    printf("\tWriting output ---- : %6.3lf\n",
           (double)(timePoint[3] - timePoint[2]) / CLOCKS_PER_SEC * 1000);
    printf("\nTotal time elapsed: %6.3lf ms\n",
           (double)(timePoint[3]) / CLOCKS_PER_SEC * 1000);

    if (argc == 4) {
        FILE *fcsv = fopen(argv[3], "a+");
        if (fcsv == NULL) {
            printf("Error opening output file.\n");
        } else {
            fprintf(
                fcsv, "%6.3lf, %6.3lf, %6.3lf, %6.3lf\n",
                (double)(timePoint[1] - timePoint[0]) / CLOCKS_PER_SEC * 1000,
                (double)(timePoint[2] - timePoint[1]) / CLOCKS_PER_SEC * 1000,
                (double)(timePoint[3] - timePoint[2]) / CLOCKS_PER_SEC * 1000,
                (double)(timePoint[3]) / CLOCKS_PER_SEC * 1000);
            fclose(fcsv);
        }
    }
    return 0;
}

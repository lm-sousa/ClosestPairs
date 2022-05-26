#include <float.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "closestpair.h"
#include "mergesort.h"
#include "types.h"

FILE *fi;
FILE *fo;

dimensionIndex_t dims = 2;
pointIndex_t numberOfPoints = 0;

double *dimensionData = NULL;
double *points = NULL;

#define dimensionMin(dimension) (dimensionData[2 * (dimension) + 0])
#define dimensionMax(dimension) (dimensionData[2 * (dimension) + 1])

#define getPointIndex(pointNumber) (dims * (pointNumber))

clock_t timePoint[11];
double timeDiff;

typedef struct {
    int processes;
    int id;
} MpiData;

MPI_Status st;

int main(int argc, char *argv[]) {
    timePoint[0] = clock();

    MpiData mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.id);

    if (mpi.id == 0) {
        timePoint[1] = clock();

        if (argc < 3) {
            printf("Usage: <executable> inputFile.txt outputFile.txt\n");
            return 1;
        }

        int numValuesRead;

        fi = fopen(argv[1], "r");
        if (fi == NULL) {
            printf("Error opening input file.\n");
            return -1;
        }

        fo = fopen(argv[2], "w");
        if (fi == NULL) {
            printf("Error opening output file.\n");
            fclose(fi);
            return -1;
        }

        timePoint[2] = clock();

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

        timePoint[3] = clock();

        dimensionData = (double *)malloc(2 * sizeof(double) * dims);
        if (dimensionData == NULL) {
            fclose(fi);
            fclose(fo);
            printf("Error allocating memory for dimension data.\n");
            return -1;
        }

        points = (double *)malloc(sizeof(double) * dims * numberOfPoints);
        if (points == NULL) {
            fclose(fi);
            fclose(fo);
            free(dimensionData);
            printf("Error allocating memory for points.\n");
            return -1;
        }

        timePoint[4] = clock();

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

        timePoint[5] = clock();

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

        timePoint[6] = clock();
    }
    MPI_Bcast(&dims, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numberOfPoints, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    int partialArraySize[mpi.processes];
    int partialArraySize2[mpi.processes];
    int partialArrayDisplacement[mpi.processes];
    for (int i = 0; i < mpi.processes; i++) {
        partialArraySize[i] = (numberOfPoints / mpi.processes +
                               (i < (numberOfPoints % mpi.processes))) *
                              dims;
        partialArraySize2[i] = partialArraySize[i];

        if (i == 0) {
            partialArrayDisplacement[i] = 0;
        } else {
            partialArrayDisplacement[i] =
                partialArrayDisplacement[i - 1] + partialArraySize[i - 1];
        }
    }

    double *partialArray =
        (double *)malloc(sizeof(double) * partialArraySize[mpi.id]);
    if (partialArray == NULL) {
        fclose(fi);
        fclose(fo);
        free(dimensionData);
        free(points);
        printf("Error allocating 'partialArray'.\n");
        return -1;
    }

    double *partialArray2 =
        (double *)malloc(sizeof(double) * partialArraySize[mpi.id]);
    if (partialArray2 == NULL) {
        fclose(fi);
        fclose(fo);
        free(dimensionData);
        free(points);
        free(partialArray);
        printf("Error allocating 'partialArray2'.\n");
        return -1;
    }

    MPI_Scatterv(points, partialArraySize, partialArrayDisplacement, MPI_DOUBLE,
                 partialArray, partialArraySize[mpi.id], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);

    mergeSort(partialArray, partialArray2, partialArraySize[mpi.id], dims, 0);
    free(partialArray2);

    for (int i = 2; i <= mpi.processes; i <<= 1) {
        int currentSize = partialArraySize[mpi.id];
        for (int j = 0; j < mpi.processes; j += i) {
            partialArraySize[j] += partialArraySize[j + (i / 2)];
        }

        if (mpi.id % i == (i / 2)) {
            // sent to id - (i/2)
            MPI_Send(partialArray, partialArraySize[mpi.id], MPI_DOUBLE,
                     mpi.id - (i / 2), 0, MPI_COMM_WORLD);
            free(partialArray);
        } else if (mpi.id % i == 0) {
            // recv from id + (i/2)
            partialArray = (double *)realloc(
                partialArray, sizeof(double) * partialArraySize[mpi.id]);
            if (partialArray == NULL) {
                fclose(fi);
                fclose(fo);
                free(dimensionData);
                free(points);
                printf("Error reallocating 'partialArray'.\n");
                return -1;
            }

            double *newArray =
                (double *)malloc(sizeof(double) * partialArraySize[mpi.id]);
            if (newArray == NULL) {
                fclose(fi);
                fclose(fo);
                free(dimensionData);
                free(points);
                free(partialArray);
                printf("Error allocating 'newArray'.\n");
                return -1;
            }

            MPI_Recv(&partialArray[currentSize],
                     partialArraySize[mpi.id + (i / 2)], MPI_DOUBLE,
                     mpi.id + (i / 2), MPI_ANY_TAG, MPI_COMM_WORLD, &st);

            pointIndex_t leftPoint = 0;
            pointIndex_t middlePoint = (currentSize) / dims - 1;
            pointIndex_t rightPoint = (partialArraySize[mpi.id]) / dims - 1;

            merge(partialArray, newArray, leftPoint, middlePoint, rightPoint,
                  dims, 0);

            free(newArray);
        }
    }

    if (mpi.id == 0) {
        free(points);
        points = partialArray;
    }

    timePoint[7] = clock();

    ////////////////////////////////
    // Compute the closest pair
    ////////////////////////////////

    partialArray = (double *)malloc(sizeof(double) * partialArraySize2[mpi.id]);
    if (partialArray == NULL) {
        fclose(fi);
        fclose(fo);
        free(dimensionData);
        free(points);
        printf("Error allocating 'partialArray'.\n");
        return -1;
    }

    MPI_Scatterv(points, partialArraySize2, partialArrayDisplacement,
                 MPI_DOUBLE, partialArray, partialArraySize2[mpi.id],
                 MPI_DOUBLE, 0, MPI_COMM_WORLD);

    BestPair bestPair;
    MinDistance minDistance;
    closestPair(partialArray, bestPair.self, &minDistance.self,
                partialArraySize2[mpi.id] / dims, dims);
    bestPair.self[0] += partialArrayDisplacement[mpi.id] / dims;
    bestPair.self[1] += partialArrayDisplacement[mpi.id] / dims;

    int iter = 0;
    printf("%d;%d: %lu-%lu --> %lf\n", mpi.id, iter++, bestPair.self[0],
           bestPair.self[1], minDistance.self);

    for (int i = 2; i <= mpi.processes; i <<= 1) {
        int currentSize = partialArraySize2[mpi.id];

        MPI_Allreduce(&minDistance.self, &minDistance.best, 1, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);

        if (mpi.id % i == (i / 2)) {
            // sent to id - (i/2)
            int sendSize = 0;
            while (sendSize < currentSize) {
                if ((partialArray[sendSize] - partialArray[getPointIndex(0)]) <
                    minDistance.best) {
                    sendSize += dims;
                } else {
                    break;
                }
            }
            MPI_Send(&sendSize, 1, MPI_INTEGER, mpi.id - (i / 2), 0,
                     MPI_COMM_WORLD);
            MPI_Send(bestPair.self, 2, MPI_LONG, mpi.id - (i / 2), 1,
                     MPI_COMM_WORLD);
            MPI_Send(&minDistance.self, 1, MPI_DOUBLE, mpi.id - (i / 2), 2,
                     MPI_COMM_WORLD);
            MPI_Send(partialArray, sendSize, MPI_DOUBLE, mpi.id - (i / 2), 3,
                     MPI_COMM_WORLD);
            free(partialArray);
        } else if (mpi.id % i == 0) {
            // recv from id + (i/2)
            int recvSize;
            MPI_Recv(&recvSize, 1, MPI_INTEGER, mpi.id + (i / 2), 0,
                     MPI_COMM_WORLD, &st);
            MPI_Recv(bestPair.recv, 2, MPI_LONG, mpi.id + (i / 2), 1,
                     MPI_COMM_WORLD, &st);
            MPI_Recv(&minDistance.recv, 1, MPI_DOUBLE, mpi.id + (i / 2), 2,
                     MPI_COMM_WORLD, &st);
            partialArraySize2[mpi.id] += recvSize;

            partialArray = (double *)realloc(
                partialArray, sizeof(double) * partialArraySize2[mpi.id]);
            if (partialArray == NULL) {
                fclose(fi);
                fclose(fo);
                free(dimensionData);
                free(points);
                printf("Error reallocating 'partialArray'.\n");
                return -1;
            }

            MPI_Recv(&partialArray[currentSize], recvSize, MPI_DOUBLE,
                     mpi.id + (i / 2), 3, MPI_COMM_WORLD, &st);

            pointIndex_t leftPointIndex =
                getPointIndex((currentSize / dims) - 1);
            pointIndex_t middlePointIndex = leftPointIndex;
            pointIndex_t rightPointIndex =
                getPointIndex((partialArraySize2[mpi.id] / dims) - 1);
            while (leftPointIndex > 0) {
                if ((partialArray[middlePointIndex] -
                     partialArray[leftPointIndex]) < minDistance.best) {
                    leftPointIndex -= dims;
                } else {
                    break;
                }
            }

            /*
            printf("%d;%d c:%d, r:%d, t:%d ---", mpi.id, iter, currentSize,
                   recvSize, partialArraySize2[mpi.id]);
            for (int i = 0; i < partialArraySize2[mpi.id]; i += dims) {
                printf(" %lE", partialArray[i]);
            }
            printf("\n");
            printf("%d;%d l:%lu, m:%lu, r:%lu\n", mpi.id, iter, leftPointIndex,
                   middlePointIndex, rightPointIndex);
            fflush(stdout);
            */

            pointIndex_t stripBestPair[2];
            double stripMinDistance;

            closestPair(&partialArray[leftPointIndex], stripBestPair,
                        &stripMinDistance,
                        (rightPointIndex - leftPointIndex) / dims + 1, dims);

            if (stripMinDistance < minDistance.self) {
                minDistance.self = stripMinDistance;
                bestPair.self[0] =
                    stripBestPair[0] +
                    (leftPointIndex + partialArrayDisplacement[mpi.id]) / dims;
                bestPair.self[1] =
                    stripBestPair[1] +
                    (leftPointIndex + partialArrayDisplacement[mpi.id]) / dims;
            }

            if (minDistance.recv < minDistance.self) {
                minDistance.self = minDistance.recv;
                bestPair.self[0] = bestPair.recv[0];
                bestPair.self[1] = bestPair.recv[1];
            }

            /*
            printf("%d;%d ---", mpi.id, iter++);
            for (int i = leftPointIndex; i <= rightPointIndex; i += dims) {
                printf(" %lE", partialArray[i]);
            }
            printf("\n");
            fflush(stdout);
            */

            printf("%d;%d: %lu-%lu --> %lf\n", mpi.id, iter++, bestPair.self[0],
                   bestPair.self[1], minDistance.self);
        }
    }

    if (mpi.id == 0) {
        free(partialArray);
    }

    timePoint[8] = clock();

    if (mpi.id == 0) {

        for (uint8_t i = 0; i < 2; i++) {
            for (dimensionIndex_t dim = 0; dim < dims; dim++) {
                fprintf(fo, "%lf ",
                        points[getPointIndex(bestPair.self[i]) + dim]);
            }
            fprintf(fo, "\n");
        }

        timePoint[9] = clock();

        fclose(fi);
        fclose(fo);
        free(dimensionData);
        free(points);

        timePoint[10] = clock();

        printf("\nD: %hu; N: %lu\n", dims, numberOfPoints);
        printf("\nTime elapsed in miliseconds:\n");
        printf("\tProgram setup ----- : %6.3lf\n",
               (double)(timePoint[0]) / CLOCKS_PER_SEC * 1000);
        printf("\tInitializing MPI -- : %6.3lf\n",
               (double)(timePoint[1] - timePoint[0]) / CLOCKS_PER_SEC * 1000);
        printf("\tOpening files ----- : %6.3lf\n",
               (double)(timePoint[2] - timePoint[1]) / CLOCKS_PER_SEC * 1000);
        printf("\tReading header ---- : %6.3lf\n",
               (double)(timePoint[3] - timePoint[2]) / CLOCKS_PER_SEC * 1000);
        printf("\tAllocating memory - : %6.3lf\n",
               (double)(timePoint[4] - timePoint[3]) / CLOCKS_PER_SEC * 1000);
        printf("\tReading dimensions  : %6.3lf\n",
               (double)(timePoint[5] - timePoint[4]) / CLOCKS_PER_SEC * 1000);
        printf("\tReading points ---- : %6.3lf\n",
               (double)(timePoint[6] - timePoint[5]) / CLOCKS_PER_SEC * 1000);
        printf("\tSorting points ---- : %6.3lf\n",
               (double)(timePoint[7] - timePoint[6]) / CLOCKS_PER_SEC * 1000);
        printf("\tFinding closest pair: %6.3lf\n",
               (double)(timePoint[8] - timePoint[7]) / CLOCKS_PER_SEC * 1000);
        printf("\tWriting output file : %6.3lf\n",
               (double)(timePoint[9] - timePoint[8]) / CLOCKS_PER_SEC * 1000);
        printf("\tCleaning up ------- : %6.3lf\n",
               (double)(timePoint[10] - timePoint[9]) / CLOCKS_PER_SEC * 1000);
        printf("\nTotal time elapsed: %6.3lf ms\n",
               (double)(timePoint[10]) / CLOCKS_PER_SEC * 1000);
    }

    MPI_Finalize();
    return 0;
}

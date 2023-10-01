#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

double computePartialSum(int start, int end) {
    double sum = 0.0;
    for (int i = start; i <= end; i++) {
        double term = 1.0 / (2.0 * i - 1);
        if (i % 2 == 0) {
            sum -= term;
        }
        else {
            sum += term;
        }
    }
    return sum;
}

double getElapsedTime(clock_t start_time, clock_t end_time) {
    double elapsedSeconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    return elapsedSeconds;
}

int main(int argc, char* argv[]) {
    int rank, size;
    int totalThreads;
    int totalOperations;
    double localSum = 0.0;
    double globalSum = 0.0;
    clock_t start_time, end_time;
    double elapsedTime;
    const double actualPi = 3.14159265358979323846;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Enter the total number of operations: ");
        fflush(stdout);
        scanf("%d", &totalOperations);

        printf("Enter the total number of threads (1 to 64): ");
        fflush(stdout);
        scanf("%d", &totalThreads);
    }

    MPI_Bcast(&totalOperations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&totalThreads, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int termsPerThread = totalOperations / totalThreads;
    int start = rank * termsPerThread + 1;
    int end = start + termsPerThread - 1;

    if (rank == totalThreads - 1) {
        end = totalOperations;
    }

    start_time = clock();
    localSum = computePartialSum(start, end);
    end_time = clock();

    MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi = 4.0 * globalSum;
        elapsedTime = getElapsedTime(start_time, end_time);
        printf("Approximate value of pi: %.15f\n", pi);
        printf("Computation time: %.15f seconds\n", elapsedTime);
        double difference = fabs(M_PI - pi);
        printf("Difference from actual pi: %.15f\n", difference);
    }

    MPI_Finalize();

    return 0;
}
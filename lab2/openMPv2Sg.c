#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define N 2200

const double epsilon = 0.00001;
const double tao = 0.0005;

void printVector(double *vector) {
	for (size_t i = 0; i < N; ++i) {
		printf("%f ", vector[i]);
	}
	printf("\n");
}

void setZeroVector(double *vector) {
	memset(vector, 0, N * sizeof(double));
}

int main(int argc, char *argv[]) {
	if (argc == 0) exit(-1);
	
	double *matrixA = malloc(N * N * sizeof(double));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			matrixA[i * N + j] = 1.0;
			if (i == j) {
				matrixA[i * N + j] = 2.0;
			}
		}
	}

	double *vectorU = calloc(sizeof(double), N);
	for (size_t i = 0; i < N; ++i) {
		vectorU[i] = sin(2 * PI * (i + 1) / N);
	}
	// printVector(vectorU);
	// printf("\n");

	double *vectorX = calloc(sizeof(double), N);
	double *vectorB = calloc(sizeof(double), N);
	double *vectorAxn_b = calloc(sizeof(double), N);

	volatile double startTime = omp_get_wtime(); 

	double normAx_b = 0, normB = 0;
	#pragma omp parallel
	{
		#pragma omp single
		setZeroVector(vectorB);
		
		#pragma omp for schedule(guided, atoi(argv[1]))
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = 0; j < N; ++j) {
				vectorB[i] += matrixA[i * N + j] * vectorU[j];
			}
		}

		#pragma omp for schedule(guided, atoi(argv[1])) reduction(+:normB)
		for (size_t i = 0; i < N; ++i) {
			double a = vectorB[i];
			normB += (a * a);
		}

		do {
			#pragma omp for schedule(guided, atoi(argv[1]))
			for (size_t i = 0; i < N; ++i) {
				vectorAxn_b[i] = 0;
				for (size_t j = 0; j < N; ++j) {
					vectorAxn_b[i] += matrixA[i * N + j] * vectorX[j];
				}
			}

			#pragma omp for schedule(guided, atoi(argv[1]))
			for (size_t i = 0; i < N; ++i) {
				vectorAxn_b[i] -= vectorB[i];
			}
			#pragma omp master
			normAx_b = 0;
			
			#pragma omp for schedule(guided, atoi(argv[1])) reduction(+:normAx_b)
			for (size_t i = 0; i < N; ++i) {
				double a = vectorAxn_b[i];
				normAx_b += (a * a);
			}

			#pragma omp for schedule(guided, atoi(argv[1])) 
			for (size_t i = 0; i < N; ++i) {
				vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
			}

		} while(normAx_b / normB > epsilon * epsilon);

	}

	volatile size_t endTime = omp_get_wtime();
    printf("%f\n", endTime - startTime);
	// printVector(vectorX);

    free(matrixA);
	free(vectorU);
	free(vectorX);
	free(vectorB);
	free(vectorAxn_b);
	return 0;
}

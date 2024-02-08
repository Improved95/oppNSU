#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define N 3

// Ax = b

void printMatrix(double *matrix) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			printf("%f ", matrix[i * N + j]);
		}
		printf("\n");
	}
}

void printVector(double *vector) {
	for (size_t i = 0; i < N; ++i) {
		printf("%f ", vector[i]);
	}
	printf("\n");
}

// for (size_t i = 0; i < N; i++) {
    //     for (size_t k = 0; k < N; k++) {
    //         for (size_t j = 0; j < N; j++) {
    //             temp[i][j] += (*this)[i][k] * source[k][j];
	// 		}
    //     }
    // }

int main() {
	static const double epsilon = 1 / 10 * 10 * 10 * 10 * 10;

	double *matrixA = malloc(N * N * sizeof(double));
	for (size_t i = 0; i < N * N; ++i) {
		matrixA[i] = 1.0;
	}
	for (size_t i = 0; i < N * N; i += N + 1) {
		matrixA[i] = 2.0;
 	}

	double *vectorU = calloc(sizeof(double), N);
	for (size_t i = 0; i < N; ++i) {
		vectorU[i] = sin((2 * PI * i) / N);
	}

	double *vectorX = calloc(sizeof(double), N);
	double *vectorB = calloc(sizeof(double), N);

	double *vectorE = calloc(sizeof(double), N); // for check less epsilon

	//take vector b
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			vectorB[j] += matrixA[i * N + j] * vectorU[j];
		}
	}

	while(1) {
		for (size_t i = 0; i < N; ++i) {
			vectorE[i] = 0;
			for (size_t j = 0; j < N; ++j) {
				vectorE[j] += matrixA[i * N + j] * vectorU[j];
			}
		}
		for (size_t i = 0; i < N; ++i) {
			vectorE[i] -= vectorB[i];
		}
		double numerator = 0, denominator = 0;
		for (size_t i = 0; i < N; ++i) {
			numerator += pow(vectorE[i], 2);
			denominator += pow(vectorB[i], 2);
		}
		numerator = sqrt(numerator);
		denominator = sqrt(denominator);
		if (numerator / denominator < epsilon) {
			break;
		}

		 

	}

	// printMatrix(matrixA);
	// printVector(vectorU);
	// printVector(vectorA);
	// printVector(vectorB);

	return 0;
}

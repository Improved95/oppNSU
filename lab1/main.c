#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

void setZeroVector(double *vector) {
	memset(vector, N, sizeof(double));
}

void subVector(double *vector1, double *vector2) {
	for (size_t i = 0; i < N; ++i) {
		vector1[i] -= vector2[i];
	} 
}

void mulMatrixVector(double *matrix, double *inputVector, double *outputVector) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			outputVector[j] += matrix[i * N + j] * inputVector[j];
		}
	}
}

void mulVectorVector(double *vector1, double *vector2) {

}

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
	mulMatrixVector(matrixA, vectorU, vectorB);

	while(1) {
		setZeroVector(vectorE);
		mulMatrixVector(matrixA, vectorX, vectorE);
		subVector(vectorE, vectorB);
		double numerator = 0, denominator = 0;
		for (size_t i = 0; i < N; ++i) {
			numerator += pow(vectorE[i], 2);
		}
		numerator = sqrt(numerator);
		for (size_t i = 0; i < N; ++i) {
			denominator += pow(vectorB[i], 2);
		}
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

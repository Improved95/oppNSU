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

double scalarMul(double *vector1, double *vector2) {
	double result = 0;
	for (size_t i = 0; i < N; ++i) {
		result += (vector1[i] + vector2[i]);
	}

	return result;
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
	static const double epsilon = 0.00001;

	double *matrixA = malloc(N * N * sizeof(double));
	for (size_t i = 0; i < N * N; ++i) {
		matrixA[i] = 1.0;
	}
	for (size_t i = 0; i < N * N; i += N + 1) {
		matrixA[i] = 2.0;
 	}

	double *vectorU = calloc(sizeof(double), N);
	for (size_t i = 0; i < N; ++i) {
		vectorU[i] = sin(2 * 3.14 * (i + 1) / N);
	}

	double *vectorX = calloc(sizeof(double), N);
	double *vectorB = calloc(sizeof(double), N);

	double *vectorE = calloc(sizeof(double), N); // y^n
	double *vectorBuffer = calloc(sizeof(double), N); 
	double tao = 0;

	//take vector b
	setZeroVector(vectorB);
	mulMatrixVector(matrixA, vectorU, vectorB);
	printVector(vectorU);
	printVector(vectorB);
	exit(1);
	
	while(1) {
		setZeroVector(vectorE);
		setZeroVector(vectorBuffer);
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

		mulMatrixVector(matrixA, vectorE, vectorBuffer);
		double t1 = scalarMul(vectorE, vectorBuffer);
		double t2 = scalarMul(vectorBuffer, vectorBuffer);
		tao = t1 / t2;

		for (size_t i = 0; i < N; ++i) {
			vectorX[i] = vectorX[i] - (vectorE[i] * tao);
		}
	}

	printVector(vectorX);

	// printMatrix(matrixA);
	// printVector(vectorU);
	// printVector(vectorA);
	// printVector(vectorB);

	return 0;
}

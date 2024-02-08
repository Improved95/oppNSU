#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M_PI 3.14159265358979323846

#define N 10
static const double epsilon = 1 / 10 * 10 * 10 * 10 * 10;

// Ax = b

void printMatrix(double *matrix) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			printf("%f ", matrix[i * N + j]);
		}
		printf("\n");
	}
}

void setMatrix(double *matrix) {
	for (size_t i = 0; i < N * N; ++i) {
		matrix[i] = 1.0;
	}
	for (size_t i = 0; i < N * N; i += N + 1) {
		matrix[i] = 2.0;
 	}
}

void setVectorU(double *vector) {
	for (size_t i = 0; i < N; ++i) {
		vector[i] = sin((2 * M_PI * i) / N);
	}
}

int charIsAchieve() {
	return 1;
}

int main() {
	double *matrix = malloc(N * N * sizeof(double));
	double *vectorA = calloc(sizeof(double), N);
	double *vectorB = calloc(sizeof(double), N);
	double *vectorU = calloc(sizeof(double), N);

	setMatrix(matrix);
	setVectorU(vectorU);
	

	while (!charIsAchieve()) {

	}

	return 0;
}

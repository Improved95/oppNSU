#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
static const double epsilon = 1 / 10 * 10 * 10 * 10 * 10;

void setMatrix(double *matrix) {
	for (size_t i = 0; i < N * N; ++i) {
		matrix[i] = 1.0;
	}
	for (size_t i = 0; i < N * N; i += N) {
		matrix[i] = 2.0;
 	}
}

int charIsAchieve() {
	return 1;
}

int main() {
	double *matrix = malloc(N * N * sizeof(double));
	double *vectorA = malloc(N * sizeof(double));
	double *vectorB = malloc(N * sizeof(double));

	while (!charIsAchieve()) {

	}

	return 0;
}

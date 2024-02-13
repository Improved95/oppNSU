#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define N 10
static int rank, sizeProc;

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

void setZeroVector(double *vector) {
	memset(vector, 0, N * sizeof(double));
}

void subVector(double *vector1, double *vector2) {
	for (size_t i = 0; i < N; ++i) {
		vector1[i] -= vector2[i];
	} 
}

double scalarMul(double *vector1, double *vector2) {
	double result = 0;
	for (size_t i = 0; i < N; ++i) {
		result += (vector1[i] * vector2[i]);
	}
	return result;
}

void mulMatrixVector(double *matrix, double *inputVector, double *outputVector) {
	for (size_t i = 0; i < N; ++i) {
			
			for (size_t j = 0; j < N; ++j) {
				outputVector[i] += matrix[i * N + j] * inputVector[j];
			}

	}
}

void mulVectorVector(double *vector1, double *vector2) {

}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const double epsilon = 0.00001;

	double *matrixA = NULL;
	if (rank == 0) {

		matrixA = сalloc(N * N, sizeof(double));
		for (size_t i = 0; i < N; ++i) {
			matrixA[i] = 1.0;
			for (size_t j = 0; i < N; ++j) {
				matrixA[i] = 2.0;
			}
		}

	}

	double *vectorU = calloc(sizeof(double), N);
	if (rank == 0) {

		for (size_t i = 0; i < N; ++i) {
			vectorU[i] = sin(2 * PI * (i + 1) / N);
		}
		printVector(vectorU);
		printf("\n");

	}

	double *vectorX = calloc(sizeof(double), N);
	double *vectorB = calloc(sizeof(double), N);

	double *vectorAxn_b = calloc(sizeof(double), N);
	double *vectorAyn = calloc(sizeof(double), N);
	double tao = 0, t1 = 0, t2 = 0;

	if (rank == 0) {

		setZeroVector(vectorB);
		mulMatrixVector(matrixA, vectorU, vectorB);

	}
	
	for(size_t k = 0; 1; ++k) {
		setZeroVector(vectorAxn_b);
		mulMatrixVector(matrixA, vectorX, vectorAxn_b);
		subVector(vectorAxn_b, vectorB);

		double numerator = 0, denominator = 0;

		for (size_t i = 0; i < N; ++i) {
			double a = vectorAxn_b[i];
			numerator += (a * a);
		}
		numerator = sqrt(numerator);

		for (size_t i = 0; i < N; ++i) {
			double a = vectorB[i];
			denominator += (a * a);
		}
		denominator = sqrt(denominator);

		if (numerator / denominator < epsilon) {
			break;
		}

		setZeroVector(vectorAyn);
		mulMatrixVector(matrixA, vectorAxn_b, vectorAyn);
		t1 = scalarMul(vectorAxn_b, vectorAyn);
		t2 = scalarMul(vectorAyn, vectorAyn);
		tao = t1 / t2;

		for (size_t i = 0; i < N; ++i) {
			vectorX[i] = vectorX[i] - (vectorAxn_b[i] * tao);
		}
	}

	printVector(vectorX);
	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define N 6
static int rank, sizeProccess;

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

void mulMatrixVector(double *pieceVector, double *inputVector, double *outputVector) {

	/*for (size_t i = 0; i < N; ++i) {
			
			for (size_t j = 0; j < N; ++j) {
				outputVector[i] += matrix[i * N + j] * inputVector[j];
			}

	}*/

	int vectorQuantity = N / sizeProccess;
	if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
		vectorQuantity++;
	}

	if (rank == 0) {
		for (size_t i = 1; i < sizeProccess; ++i) {
			MPI_Send(inputVector, N, MPI_DOUBLE, i, 1991, MPI_COMM_WORLD);
		}
	}
	if (rank != 0) {
			MPI_Recv();
	}
}

void mulVectorVector(double *vector1, double *vector2) {

}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const double epsilon = 0.00001;

	double *vectorPieceOfMatrix = NULL;
	int vectorQuantity = N / sizeProccess;
	if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
		vectorQuantity++;
	}

	vectorPieceOfMatrix = calloc(vectorQuantity * N, sizeof(double));
	for (size_t i = 0; i < vectorQuantity; ++i) {
		for (size_t j = 0; j < N; ++j) {
			vectorPieceOfMatrix[i * N + j] = 1;
			if ((i * N + j) % N == i + rank + 1 * rank) {
				vectorPieceOfMatrix[i * N + j] = 2;
			}
		}
	}

	double *vectorU = NULL;
	if (rank == 0) {

		vectorU = calloc(sizeof(double), N);
		for (size_t i = 0; i < N; ++i) {
			vectorU[i] = sin(2 * PI * (i + 1) / N);
		}
		// printVector(vectorU);
		// printf("\n");

	}

	double *vectorX = NULL;
	double *vectorB = NULL;
	if (rank == 0) {

		vectorX = calloc(sizeof(double), N);
		vectorB = calloc(sizeof(double), N);

	}

	double *vectorAxn_b = NULL;
	double *vectorAyn = NULL;
	if (rank == 0) {

		vectorAxn_b = calloc(sizeof(double), N);
		vectorAyn = calloc(sizeof(double), N);

	}
	double tao = 0.00025;

	if (rank == 0) {

		setZeroVector(vectorB);
		mulMatrixVector(vectorPieceOfMatrix, vectorU, vectorB);

	}
	
	int isComplete = 0;
	for(size_t k = 0; 1; ++k) {
		if (rank == 0) {

			setZeroVector(vectorAxn_b);

		}

		// mulMatrixVector(matrixA, vectorX, vectorAxn_b, vectorPieceOfMatrix);

		if (rank == 0) {

			subVector(vectorAxn_b, vectorB);

		}
			

		if (rank == 0) {

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
				isComplete = 1;
			}
			
			for (size_t i = 1; i < sizeProccess; ++i) {
				MPI_Send(&isComplete, 1, MPI_INT, i, 199, MPI_COMM_WORLD);
			}
		}

		if (rank != 0) {
			MPI_Recv(&isComplete, 1, MPI_INT, 0, 199, MPI_COMM_WORLD, &st);
		}

		if (isComplete) break;


		//--------------------------==---------------------------

		for (size_t i = 0; i < N; ++i) {
			vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
		}
	}

	printVector(vectorX);

	free(vectorPieceOfMatrix);
	free(vectorU);
	free(vectorX);
	free(vectorB);
	free(vectorAxn_b);
	return 0;
}
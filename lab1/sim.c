#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define N 500

static int rank, sizeProccess;
const double epsilon = 0.00001;
const double tao = 0.0003;

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

void mulMatrixVector(double *pieceVector, double *inputVector, double *outputVector, 
						double *vectorBuffer, MPI_Status st) {

	int vectorQuantity = N / sizeProccess;
	if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
		vectorQuantity++;
	}

	if (rank != 0) {

		inputVector = vectorBuffer;
		
	}
	MPI_Bcast(inputVector, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {

		for (size_t i = 0; i < vectorQuantity; ++i) {
			for (size_t j = 0; j < N; ++j) {
				outputVector[i] += pieceVector[i * N + j] * inputVector[j];
			}
		}

	}
	if (rank != 0) {

		double res = 0;
		for (size_t i = 0; i < vectorQuantity; ++i) {		
			res = 0;												
			for (size_t j = 0; j < N; ++j) {
				res += pieceVector[i * N + j] * inputVector[j];
			}

			MPI_Send(&res, 1, MPI_DOUBLE, 0, 1992, MPI_COMM_WORLD);
		}

	}

	if (rank == 0) {

		double res = 0;
		size_t posInOutputVector = vectorQuantity;
		for (size_t i = 1; i < sizeProccess; ++i) {
			int vectorQuantityInAnotherProcess = N / sizeProccess;
			if ((N % sizeProccess != 0) && (N % sizeProccess >= i + 1)) {
				vectorQuantityInAnotherProcess++;
			}
			
			for (size_t j = 0; j < vectorQuantityInAnotherProcess; ++j) {
				MPI_Recv(&res, 1, MPI_DOUBLE, i, 1992, MPI_COMM_WORLD, &st);
				outputVector[posInOutputVector] = res;
				++posInOutputVector;
			}
		}
		
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *pieceVector = NULL;
	int vectorQuantity = N / sizeProccess;
	if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
		vectorQuantity++;
	}
	int vectorQuantityInPrevProcess = 0;
	for (size_t i = 0; i < rank; ++i) {
		vectorQuantityInPrevProcess += N / sizeProccess;
		if ((N % sizeProccess != 0) && (N % sizeProccess >= i + 1)) {
			vectorQuantityInPrevProcess++;
		}
	}
	pieceVector = calloc(vectorQuantity * N, sizeof(double));
	for (size_t i = 0; i < vectorQuantity; ++i) {
		for (size_t j = 0; j < N; ++j) {
			pieceVector[i * N + j] = 1;
			if (j == i + vectorQuantityInPrevProcess) {
				pieceVector[i * N + j] = 2;
			}
		}
	}

	double *vectorU = NULL;
	if (rank == 0) {

		printf("MPI\n");
		vectorU = calloc(N, sizeof(double));
		for (size_t i = 0; i < N; ++i) {
			vectorU[i] = sin(2 * PI * (i + 1) / N);
		}
		// printVector(vectorU);
		// printf("\n");

	}

	double *vectorX = NULL;
	double *vectorB = NULL;
	if (rank == 0) {

		vectorX = calloc(N, sizeof(double));
		vectorB = calloc(N, sizeof(double));

	}

	double *vectorBuffer = NULL;
	if (rank != 0) {

		vectorBuffer = calloc(N, sizeof(double));

	}

	double *vectorAxn_b = NULL;
	if (rank == 0) {

		vectorAxn_b = calloc(N, sizeof(double));

	}

	if (rank == 0) {

		setZeroVector(vectorB);

	}
	
	mulMatrixVector(pieceVector, vectorU, vectorB, vectorBuffer, st);

	double startTime = 0;
	if (rank == 0) {

		startTime = MPI_Wtime();

	}

	int isComplete = 0;
	for(size_t k = 0; 1; ++k) {
		if (rank == 0) {

			setZeroVector(vectorAxn_b);

		}

		mulMatrixVector(pieceVector, vectorX, vectorAxn_b, vectorBuffer, st);

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

		if (isComplete) {
			if (rank == 0) {

				printf("%ld\n", k);
			
			}
			break; 
		}

		if (rank == 0) {

			for (size_t i = 0; i < N; ++i) {
				vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
			}
		
		}
	}

	if (rank == 0) {
		
		double endTime = MPI_Wtime();
		printf("%f\n", endTime - startTime);
		// printVector(vectorX);
	
	}
	free(pieceVector);
	free(vectorBuffer);
	free(vectorU);
	free(vectorX);
	free(vectorB);
	free(vectorAxn_b);

	MPI_Finalize();
	return 0;
}

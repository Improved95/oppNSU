#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define N 2200

static int rank, sizeProccess;
const double epsilon = 0.00001;
const double tao = 0.0005;

void printMatrix(double *matrix) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			printf("%f ", matrix[i * N + j]);
		}
		printf("\n");
	}
}

void breakProgramm() {
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(-1);
}

void printVector(double *vector) {
	for (size_t i = 0; i < N; ++i) {
		printf("%f ", vector[i]);
	}
	printf("\n");
}

void printVectorv2(double *vector, size_t vectorSize) {
	for (size_t i = 0; i < (size_t)sizeProccess; ++i) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (i == (size_t)rank) {
			// printf("rank %d: ", rank);
			for (size_t j = 0; j < vectorSize; ++j) {
				printf("%f ", vector[j]);
			}
			// printf("\n");
		}
	}
	if (rank == sizeProccess - 1) {
			printf("\n");
			printf("\n");
	}
}

void setZeroVector(double *vector) {
	memset(vector, 0, N * sizeof(double));
}

void setZeroVectorV2(double *vector, size_t vectorSize) {
	memset(vector, 0, vectorSize * sizeof(double));
}

void subVector(double *vector1, double *vector2) {
	for (size_t i = 0; i < N; ++i) {
		vector1[i] -= vector2[i];
	} 
}

void subVectorV2(double *vector1, double *vector2, size_t sizeVector) {
	for (size_t i = 0; i < sizeVector; ++i) {
		vector1[i] -= vector2[i];
	} 
}

double getNorm(double *vector, size_t sizeVector) {
	double sum = 0;
	for (size_t i = 0; i < sizeVector; ++i) {
		double a = vector[i];
		sum += (a * a);
	}

	double res = 0;
	MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return res;
}

void mulMatrixVector(double *pieceVector, double *inputVector, double *outputVector, 
		size_t vectorSizeInCurrentProcess) {

	MPI_Bcast(inputVector, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
		for (size_t j = 0; j < N; ++j) {
			outputVector[i] += pieceVector[i * N + j] * inputVector[j];
		}
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *pieceVector = NULL;
	size_t vectorSizeInCurrentProcess = N / sizeProccess;
	if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
		vectorSizeInCurrentProcess++;
	}

	size_t sumSizeVectorInPrevProcesses = 0;
	for (size_t i = 0; i < (size_t)rank; ++i) {
		sumSizeVectorInPrevProcesses += N / sizeProccess;
		if ((N % sizeProccess != 0) && (N % sizeProccess >= i + 1)) {
			sumSizeVectorInPrevProcesses++;
		}
	}

	pieceVector = calloc(vectorSizeInCurrentProcess * N, sizeof(double));
	for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
		for (size_t j = 0; j < N; ++j) {
			pieceVector[i * N + j] = 1;
			if (j == i + sumSizeVectorInPrevProcesses) {
				pieceVector[i * N + j] = 2;
			}
		}
	}

	double *vectorU = NULL;
	if (rank == 0) {

		printf("MPIv1\n");
		vectorU = calloc(N, sizeof(double));
		for (size_t i = 0; i < N; ++i) {
			vectorU[i] = sin(2 * PI * (i + 1) / N);
		}
		// printVector(vectorU);
		// printf("\n");

	}

	double *vectorX = NULL;
	double *completeVectorAxn_b = NULL;
	if (rank == 0) {
		
		vectorX = calloc(N, sizeof(double));
		completeVectorAxn_b = calloc(N, sizeof(double));

	}

	double *vectorBuffer = NULL;
	if (rank != 0) {

		vectorBuffer = calloc(N, sizeof(double));

	}

	double *vectorB = calloc(vectorSizeInCurrentProcess, sizeof(double));
	double *vectorAxn_b = calloc(vectorSizeInCurrentProcess, sizeof(double));

	if (rank == 0) {

		setZeroVector(vectorB);

	}
	
	if (rank != 0) {

		vectorU = vectorBuffer;

	}

	double startTime = MPI_Wtime();

	mulMatrixVector(pieceVector, vectorU, vectorB, vectorSizeInCurrentProcess);

	double normB = getNorm(vectorB, vectorSizeInCurrentProcess);

	if (rank != 0) {

			vectorX = vectorBuffer;
		
	}

	int *recvCounts = calloc(sizeProccess, sizeof(int));
	int *displs = calloc(sizeProccess, sizeof(int));
	size_t displsCounter = 0;
	for (size_t i = 0; i < (size_t)sizeProccess - 1; ++i) {
		displsCounter += N / sizeProccess;
		recvCounts[i] = N / sizeProccess;
		if ((N % sizeProccess != 0) && (N % sizeProccess >= i + 1)) {
			displsCounter++;
			recvCounts[i]++;
		}
		displs[i + 1] = displsCounter;
	}
	recvCounts[sizeProccess - 1] = N / sizeProccess;


	for(size_t k = 0; 1; ++k) {
		setZeroVectorV2(vectorAxn_b, vectorSizeInCurrentProcess);
		mulMatrixVector(pieceVector, vectorX, vectorAxn_b, vectorSizeInCurrentProcess);
		subVectorV2(vectorAxn_b, vectorB, vectorSizeInCurrentProcess);

		double normAx_b = getNorm(vectorAxn_b, vectorSizeInCurrentProcess);
		if (normAx_b / normB < epsilon * epsilon) {
			// if (rank == 0) printf("iterations: %ld\n", k);
			break;
		}

		MPI_Gatherv(vectorAxn_b, vectorSizeInCurrentProcess, MPI_DOUBLE, 
						completeVectorAxn_b, recvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		if (rank == 0) {
			
			for (size_t i = 0; i < N; ++i) {
				vectorX[i] = vectorX[i] - (tao * completeVectorAxn_b[i]);
			}
		
		}
	}

	double endTime = MPI_Wtime();
	double time = endTime - startTime;

	double finalTime = 0;
	MPI_Reduce(&time, &finalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0) {

		// printVector(vectorX);
		printf("time: %f\n", finalTime);

	}

	free(pieceVector);
	free(vectorBuffer);
	if (rank == 0) {
		
		free(vectorU);
		free(vectorX);

	}
	free(vectorB);
	free(vectorAxn_b);
	free(completeVectorAxn_b);

	MPI_Finalize();
	return 0;
}

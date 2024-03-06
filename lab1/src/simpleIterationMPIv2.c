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

void printVector(double *vector, size_t vectorSize) {
	for (size_t i = 0; i < vectorSize; ++i) {
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

void setZeroVector(double *vector, size_t vectorSize) {
	memset(vector, 0, vectorSize * sizeof(double));
}

void subVector(double *vector1, double *vector2, size_t sizeVector) {
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

void mulMatrixVector(double *pieceVector, double *circleVector, double *outputVector, double *circleBuffer,
						size_t vectorSizeInCurrentProcess, size_t shiftSize, size_t sumSizeVectorInPrevProcesses) {
	
	MPI_Status st;

	size_t indexInVector = sumSizeVectorInPrevProcesses; //позиция в строке в которой работаем
	size_t numberBlockInVector = rank; // номер блока в котором работаем
	for (size_t i = 0; i < (size_t)sizeProccess; ++i) {
		indexInVector %= N;

		for (size_t j = 0; j < vectorSizeInCurrentProcess; ++j) {
			for (size_t k = 0; k < shiftSize; ++k) {
				outputVector[j] += pieceVector[(j * N + k + indexInVector)] * circleVector[k];
			}
		}

		size_t blockSize = N / sizeProccess;
		if ((N % sizeProccess != 0) && (N % sizeProccess >= numberBlockInVector + 1)) {
			blockSize++;
		}
		
		if ((N % sizeProccess != 0) && (shiftSize > blockSize)) {
			indexInVector--;
		}
		indexInVector += shiftSize;
		numberBlockInVector++;
		numberBlockInVector %= sizeProccess;

		/*сдвиг циклический*/
		if (sizeProccess > 1) {
			if (rank % 2 == 0) {
				MPI_Send(circleVector, shiftSize, MPI_DOUBLE, (rank - 1 + sizeProccess) % sizeProccess, 121212, MPI_COMM_WORLD);
			}
			if (rank % 2 != 0) {
				memmove(circleBuffer, circleVector, shiftSize * sizeof(double));
				MPI_Recv(circleVector, shiftSize, MPI_DOUBLE, (rank + 1) % sizeProccess, 121212, MPI_COMM_WORLD, &st);
			}
			
			if (rank % 2 != 0) {
				MPI_Send(circleBuffer, shiftSize, MPI_DOUBLE, (rank - 1 + sizeProccess) % sizeProccess, 121212, MPI_COMM_WORLD);	
			}
			if (rank % 2 == 0) {
				MPI_Recv(circleVector, shiftSize, MPI_DOUBLE, (rank + 1) % sizeProccess, 121212, MPI_COMM_WORLD, &st);
			}
		}
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *pieceVectorOfMatrix = NULL;
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
	pieceVectorOfMatrix = calloc(vectorSizeInCurrentProcess * N, sizeof(double));
	for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
		for (size_t j = 0; j < N; ++j) {
			pieceVectorOfMatrix[i * N + j] = 1;
			if (j == i + sumSizeVectorInPrevProcesses) {
				pieceVectorOfMatrix[i * N + j] = 2;
			}
		}
	}

	size_t shiftSize = N / sizeProccess + (N % sizeProccess != 0);

	double *vectorU = calloc(shiftSize, sizeof(double));
	for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
		vectorU[i] = sin(2 * PI * (i + 1 + sumSizeVectorInPrevProcesses) / N);
	}

	if (rank == 0) {

		printf("MPIv2\n");

	}

	double *vectorX = calloc(shiftSize, sizeof(double));
	double *vectorB = calloc(shiftSize, sizeof(double));
	double *vectorAxn_b = calloc(shiftSize, sizeof(double));
	double *circleBuffer = NULL;
	if (rank % 2 != 0) {

		circleBuffer = calloc(shiftSize, sizeof(double));
		
	}

	double startTime = MPI_Wtime();

	mulMatrixVector(pieceVectorOfMatrix, vectorU, vectorB, circleBuffer,
						vectorSizeInCurrentProcess, shiftSize, sumSizeVectorInPrevProcesses);


	double normB = getNorm(vectorB, vectorSizeInCurrentProcess);

	for(size_t k = 0; 1; ++k) {
		setZeroVector(vectorAxn_b, vectorSizeInCurrentProcess);
		mulMatrixVector(pieceVectorOfMatrix, vectorX, vectorAxn_b, circleBuffer,
							vectorSizeInCurrentProcess, shiftSize, sumSizeVectorInPrevProcesses);
		subVector(vectorAxn_b, vectorB, vectorSizeInCurrentProcess);
			
		double normAx_b = getNorm(vectorAxn_b, vectorSizeInCurrentProcess);
		if (normAx_b / normB < epsilon * epsilon) {
			if (rank == 0) printf("iterations: %ld\n", k);
			break;
		}
		
		for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
			vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
		}
	}

	double endTime = MPI_Wtime();
	double time = endTime - startTime;

	double finalTime = 0;
	MPI_Reduce(&time, &finalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	// printVectorv2(vectorU, vectorSizeInCurrentProcess);
	// printVectorv2(vectorX, vectorSizeInCurrentProcess);

	if (rank == 0) {

		printf("%f\n", finalTime);

	}

	free(pieceVectorOfMatrix);
	free(vectorU);
	free(vectorX);
	free(vectorB);
	free(vectorAxn_b);

	MPI_Finalize();
	return 0;
}

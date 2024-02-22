#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define N 6

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

void breakProgramm() {
	MPI_Barrier(MPI_COMM_WORLD);
	exit(-1);
}

void printVector(double *vector) {
	for (size_t i = 0; i < N; ++i) {
		printf("%f ", vector[i]);
	}
	printf("\n");
}

void printVectorv2(double *vector, size_t vectorSize) {
	for (size_t i = 0; i < sizeProccess; ++i) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (i == rank) {
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
		size_t vectorSizeInCurrentProcess, size_t sumSizeVectorInPrevProcesses ) {

	MPI_Bcast(inputVector, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
		for (size_t j = 0; j < N; ++j) {
			outputVector[i] += pieceVector[i * N + j] * inputVector[j];
		}
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *pieceVector = NULL;
	size_t vectorSizeInCurrentProcess = N / sizeProccess;
	if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
		vectorSizeInCurrentProcess++;
	}
	size_t sumSizeVectorInPrevProcesses = 0;
	for (size_t i = 0; i < rank; ++i) {
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

		printf("MPI\n");
		vectorU = calloc(N, sizeof(double));
		for (size_t i = 0; i < N; ++i) {
			// vectorU[i] = sin(2 * PI * (i + 1) / N);
			vectorU[i] = i + 1;
		}
		printVector(vectorU);
		printf("\n");

	}

	double *vectorX = NULL;
	if (rank == 0) {

		vectorX = calloc(N, sizeof(double));

	}

	double *vectorBuffer = NULL;
	if (rank != 0) {

		vectorBuffer = calloc(N, sizeof(double));

	}

	double *vectorB = NULL;
	double *vectorAxn_b = NULL;
	vectorB = calloc(vectorSizeInCurrentProcess, sizeof(double));
	vectorAxn_b = calloc(vectorSizeInCurrentProcess, sizeof(double));

	if (rank == 0) {

		setZeroVector(vectorB);

	}
	
	if (rank != 0) {

		vectorU = vectorBuffer;

	}
	mulMatrixVector(pieceVector, vectorU, vectorB, 
						vectorSizeInCurrentProcess, sumSizeVectorInPrevProcesses);

	double normB = getNorm(vectorB, vectorSizeInCurrentProcess);

	double startTime = MPI_Wtime();

	if (rank != 0) {

			vectorX = vectorBuffer;
		
	}
	int isComplete = 0;
	for(size_t k = 0; 1; ++k) {
		setZeroVectorV2(vectorAxn_b, vectorSizeInCurrentProcess);
		
		mulMatrixVector(pieceVector, vectorX, vectorAxn_b,
							vectorSizeInCurrentProcess, sumSizeVectorInPrevProcesses);

		if (rank == 0) {

			subVector(vectorAxn_b, vectorB);

		}
			
		if (rank == 0) {

			double numerator = 0;
			
			for (size_t i = 0; i < N; ++i) {
				double a = vectorAxn_b[i];
				numerator += (a * a);
			}
			numerator = sqrt(numerator);

			if (numerator / normB < epsilon) {
				isComplete = 1;
			}

		}

		MPI_Bcast(&isComplete, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (isComplete) {
			break; 
		}

		if (rank == 0) {

			for (size_t i = 0; i < N; ++i) {
				vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
			}
		
		}
	}

	double endTime = MPI_Wtime();
	double time = endTime - startTime;

	double finalTime = 0;
	MPI_Reduce(&time, &finalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		
		printVector(vectorX);
		printf("%f\n", finalTime);

	}

	free(pieceVector);
	free(vectorBuffer);
	if (rank == 0) {
		
		free(vectorU);
		free(vectorX);

	}
	free(vectorB);
	free(vectorAxn_b);

	MPI_Finalize();
	return 0;
}

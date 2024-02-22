#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define N 5

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

void printVector(double *vector, size_t vectorSize) {
	for (size_t i = 0; i < vectorSize; ++i) {
		printf("%f ", vector[i]);
	}
	printf("\n");
}

void printVectorv2(double *vector, size_t vectorSize) {
	for (size_t i = 0; i < sizeProccess; ++i) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (i == rank) {
			printf("rank %d: ", rank);
			for (size_t j = 0; j < vectorSize; ++j) {
				printf("%f ", vector[j]);
			}
			printf("\n");
		}
	}
	if (rank == sizeProccess - 1) {
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

void mulMatrixVector(double *pieceVector, double *circleVector, double *outputVector,
						size_t vectorSizeInCurrentProcess, size_t shiftSize, size_t sumSizeVectorInPrevProcesses) {
	
	MPI_Request req[2];
	MPI_Status st;
	
	// for (int i = 0; i < sizeProccess; ++i) {
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// 	if (i == rank) {
	// 		printf("rank %d, vs %ld, sum %ld, shift %ld\n", rank, vectorSizeInCurrentProcess, sumSizeVectorInPrevProcesses, shiftSize);
	// 	}
	// }
	// breakProgramm();

	size_t indexInVector = sumSizeVectorInPrevProcesses; //позиция в строке в которой работаем
	size_t numberBlockInVector = rank; // номер блока в котором работаем
	for (size_t i = 0; i < sizeProccess; ++i) {
		indexInVector %= N;

		for (size_t j = 0; j < vectorSizeInCurrentProcess; ++j) {
			for (size_t k = 0; k < shiftSize; ++k) {
				outputVector[j] += pieceVector[(j * N + k + indexInVector) % (N * vectorSizeInCurrentProcess)] * circleVector[k];
				// if (rank == 0 && j == 1) {
				// 	printf("%d %d %d %ld\n", (int)outputVector[j], 
				// 	(int)pieceVector[(j * N + k + indexInVector) % (N * vectorSizeInCurrentProcess)], (int)circleVector[k], indexInVector);
				// }
			}
		}

		size_t blockSize = N / sizeProccess;
		if ((N % sizeProccess != 0) && (N % sizeProccess >= numberBlockInVector + 1)) {
			blockSize++;
		}
		// if (rank == 0) {
		// 	printf("%ld ", blockSize);
		// }
		if ((N % sizeProccess != 0) && (shiftSize > blockSize)) {
			indexInVector--;
		}
		indexInVector += shiftSize;
		numberBlockInVector++;
		numberBlockInVector %= sizeProccess;
		
		// if (rank == 0) {
			// printVector(outputVector, vectorSizeInCurrentProcess);
		// }

		/*сдвиг циклический*/
		if (sizeProccess > 1) {		   
			MPI_Isend(circleVector, shiftSize, MPI_DOUBLE, (rank - 1 + sizeProccess) % sizeProccess, 12345, MPI_COMM_WORLD, &req[0]);
			MPI_Irecv(circleVector, shiftSize, MPI_DOUBLE, (rank + 1	  	 	   ) % sizeProccess, 12345, MPI_COMM_WORLD, &req[1]);
			MPI_Waitall(2, req, &st);
		}
	}
	
	// printVectorv2(outputVector, vectorSizeInCurrentProcess);
	// breakProgramm();

}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *pieceVectorOfMatrix = NULL;
	int vectorSizeInCurrentProcess = N / sizeProccess;
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
	vectorU = calloc(N, sizeof(double));
	for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
		// vectorU[i] = sin(2 * PI * (i + 1 + sumSizeVectorInPrevProcesses) / N);
		vectorU[i] = i + sumSizeVectorInPrevProcesses + 1;
	}

	// if (rank == 0) { printf("MPIv2\n"); }

	// printVectorv2(vectorU, vectorSizeInCurrentProcess);
	// breakProgramm();

	double *vectorX = NULL;
	double *vectorB = NULL;
	// double *vectorBuffer = NULL;
	vectorX = calloc(shiftSize, sizeof(double));
	vectorB = calloc(shiftSize, sizeof(double));

	double *vectorAxn_b = NULL;
	vectorAxn_b = calloc(shiftSize, sizeof(double));

	mulMatrixVector(pieceVectorOfMatrix, vectorU, vectorB, 
						vectorSizeInCurrentProcess, shiftSize, sumSizeVectorInPrevProcesses);

	printVectorv2(vectorB, vectorSizeInCurrentProcess);
	// breakProgramm();

	double normB = getNorm(vectorB, vectorSizeInCurrentProcess);
	if (rank == 0) {
		printf("%f", normB);
	}
	breakProgramm();

	int isComplete = 0;
	for(size_t k = 0; 1; ++k) {
		setZeroVector(vectorAxn_b, vectorSizeInCurrentProcess);
		mulMatrixVector(pieceVectorOfMatrix, vectorX, vectorAxn_b, 
							vectorSizeInCurrentProcess, shiftSize, sumSizeVectorInPrevProcesses);
		subVector(vectorAxn_b, vectorB, vectorSizeInCurrentProcess);
			


		

			// for (size_t i = 0; i < N; ++i) {
			// 	vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
			// }
	}

	free(pieceVectorOfMatrix);
	free(vectorU);
	free(vectorX);
	free(vectorB);
	free(vectorAxn_b);

	MPI_Finalize();
	return 0;
}
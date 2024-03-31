#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

static int rank, sizeProccess;
const int N1 = 32;
const int N2 = 16;
const int M1 = N2;
const int M2 = N1;

void printEnter() {
    printf("\n");
}

void printMatrix(double *matrix, int n1, int n2) {
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			printf("%f ", matrix[i * n2 + j]);
		}
		printf("\n");
	}
}

void breakProgramm() {
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(-1);
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *matrix1 = NULL;
    double *matrix2 = NULL;
    double *matrix3 = NULL;
    if (rank == 0) {

        double *matrix1 = malloc(N1 * N2 * sizeof(double));
        for (int i = 0; i < N1 * N2; ++i) {
            matrix1[i] = i;
        }

        double *matrix2 = malloc(M1 * M2 * sizeof(double));
        for (int i = 0; i < N1 * N2; ++i) {
            matrix2[i] = i;
        }

        double *matrix2 = calloc(N2 * M1, sizeof(double));
    }



    MPI_Finalize();
	return 0;
}

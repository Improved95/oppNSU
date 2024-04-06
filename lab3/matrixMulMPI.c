#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>

#define DIMS_COUNT 2
#define X 0
#define Y 1

#define N1 32
#define N2 16

#define M1 24
#define M2 N1

#define DIMS_X 2
#define DIMS_Y 2

static int rank, sizeProccess;

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

void initCommunicators(const int dims[DIMS_COUNT], MPI_Comm *commGrid, MPI_Comm *commRows, MPI_Comm *commColumns) {
    int reorder = 1;
    int periods[DIMS_COUNT] = {};
    int subDims[DIMS_COUNT] = {};

    MPI_Cart_create(MPI_COMM_WORLD, DIMS_COUNT, dims, periods, reorder, commGrid);

    subDims[X] = false;
    subDims[Y] = true;
    MPI_Cart_sub(*commGrid, subDims, commRows);

    subDims[X] = true;
    subDims[Y] = false;
    MPI_Cart_sub(*commGrid, subDims, commColumns);    
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


    int dims[DIMS_COUNT] = {DIMS_X, DIMS_Y};
    MPI_Dims_create(sizeProccess, DIMS_COUNT, dims);

    MPI_Comm commGrid;
    MPI_Comm commRows;
    MPI_Comm commColumns;
    initCommunicators(dims, commGrid, commRows, commColumns);
    
    int coords[DIMS_COUNT] = {};
    MPI_Cart_coords(commGrid, rank, DIMS_COUNT, coords);

    double startTime = MPI_Wtime();

    double *matrix1Block = NULL;
    double *matrix2Block = NULL;
    double *matrix3Block = NULL;
    


    double finishTime = MPI_Wtime();

    if (rank == 0) {

        printf("Time: %f\n", finishTime - startTime);

        free(matrix1);
        free(matrix2);
        free(matrix3);

    }
    
    free(matrix1Block);
    free(matrix2Block);
    free(matrix3Block);

    MPI_Finalize();
	return 0;
}

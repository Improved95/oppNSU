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
#define N2 24 
#define N3 16

#define DIMS_X 2
#define DIMS_Y 2

static int rank, sizeProccess;

void breakProgramm() {
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(-1);
}

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

/*----------*/

void fillMatrix(double *matrix, int n1, int n2) {
    for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			matrix[i * n2 + j] = i + j;
		}
	}
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

void splitA(double *matrix1, double *matrix1Block, int matrix1BlockSize, 
            int n2, int coordsY, MPI_Comm commRows, MPI_Comm commColumns) {
    
    if (coordsY == 0) {

        MPI_Scatter(matrix1, matrix1BlockSize * n2, MPI_DOUBLE, matrix1Block, matrix1BlockSize * n2, MPI_DOUBLE, 0, commColumns);
    
    }

    MPI_Bcast(matrix1Block, matrix1BlockSize * n2, MPI_DOUBLE, 0, commRows);
}
 
int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int dims[DIMS_COUNT] = {DIMS_X, DIMS_Y};
    MPI_Dims_create(sizeProccess, DIMS_COUNT, dims);

    MPI_Comm commGrid;
    MPI_Comm commRows;
    MPI_Comm commColumns;
    initCommunicators(dims, commGrid, commRows, commColumns);
    
    int coords[DIMS_COUNT] = {};
    MPI_Cart_coords(commGrid, rank, DIMS_COUNT, coords);

    double *matrix1 = NULL;
    double *matrix2 = NULL;
    double *matrix3 = NULL;
    
    int matrix1BlockSize = ceil((double) N1 / dims[X]);
    int matrix2BlockSize = ceil((double) N3 / dims[Y]);
    int alignedN1 = matrix1BlockSize * dims[X];
    int alignedN3 = matrix2BlockSize * dims[Y];
    if (rank == 0) {

        matrix1 = malloc(alignedN1 * N2 * sizeof(double));
        matrix2 = malloc(N2 * alignedN3 * sizeof(double));
        matrix3 = malloc(alignedN1 * alignedN3 * sizeof(double));
        
        fillMatrix(matrix1, alignedN1, N2);
        fillMatrix(matrix2, N2, alignedN3);
    }

    double startTime = MPI_Wtime();

    double *matrix1Block = malloc(matrix1BlockSize * N2 * sizeof(double));
    double *matrix2Block = malloc(matrix2BlockSize * N2 * sizeof(double));
    double *matrix3Block = malloc(matrix1BlockSize * matrix2BlockSize * sizeof(double));

    splitA(matrix1, matrix1Block, matrix1BlockSize, N2, coords[Y], commRows, commColumns);
    splitB();

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
    MPI_Comm_free(&commGrid);
    MPI_Comm_free(&commRows);
    MPI_Comm_free(&commColumns);

    MPI_Finalize();
	return 0;
}

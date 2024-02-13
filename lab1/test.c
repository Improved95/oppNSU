#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

static int x;

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Status st;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size >= 2) {
        if (rank == 0) {
            x = 10;
            // MPI_Send(&x, 1, MPI_INT, 1, 01, MPI_COMM_WORLD);
        }
        if (rank == 1) {
            x = 5;
            // MPI_Recv(&x, 1, MPI_INT, 0, 01, MPI_COMM_WORLD, &st);
        }
    }

    MPI_Finalize();

    printf("Process: %d, size: %d, x: %d\n", rank, size, x);

    return 0;
}

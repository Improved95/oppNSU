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

    int a, b;
    MPI_Get_version(&a, &b);


    MPI_Finalize();

    printf("%d %d", a, b);

    return 0;
}

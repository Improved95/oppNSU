rm matrixMulMPI.out
clear
mpicc matrixMulMPI.c -lm -Wall -Werror -Wextra -O2 -o matrixMulMPI.out
mpiexec -n $1 ./matrixMulMPI.out
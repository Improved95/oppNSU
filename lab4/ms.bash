rm matrixMulMPI.out
clear
mpicc solveEquation.c -lm -Wall -Werror -Wextra -O2 -o solveEquation.out
mpiexec -n $1 ./solveEquation.out
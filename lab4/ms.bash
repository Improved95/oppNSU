rm matrixMulMPI.out
clear
mpicc solveEquation.c -lm -Wall -Werror -Wextra -O2 -o solveEquation.out

# mpiexec -n 1 ./solveEquation.out

# echo -e "\n"
# mpiexec -n 2 ./solveEquation.out

# echo -e "\n"
# mpiexec -n 4 ./solveEquation.out

# echo -e "\n"
# mpiexec -n 8 ./solveEquation.out

# echo -e "\n"
# mpiexec -n 15 ./solveEquation.out
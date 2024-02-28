mpicc src/simpleIterationMPIv$1.c -O3 -Wall -Werror -Wextra -lm -o simv$1.out
mpiexec -n $2 ./simv1.out


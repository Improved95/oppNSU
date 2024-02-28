mpicc src/simpleIterationMPIv$1.c -O3 -lm -o simv$1.out
mpiexec -n $2 ./simv$1.out


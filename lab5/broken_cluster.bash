# cmake ../../oppNSU/lab5 make && mpiexec -np 8 ./cluster
# cmake ../ && make && mpiexec -np 8 ./cluster
# make  && mpiexec -np 8 ./cluster

rm broken_cluster.out
clear

mpicc broken_cluster.c -lm -Wall -Werror -Wextra -O2 -o broken_cluster.out
mpiexec -n 8 ./broken_cluster.out
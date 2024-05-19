# cmake ../../oppNSU/lab5 make && mpiexec -np 8 ./cluster
# cmake ../ && make && mpiexec -np 8 ./cluster
# make  && mpiexec -np 8 ./cluster

rm cluster.out
clear

mpicc cluster.c -lm -Wall -Werror -Wextra -O2 -o cluster.out
mpiexec -n 8 ./cluster.out